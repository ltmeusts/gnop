/*
Copyright (c) 2025, This project is built by T. Li, Suzhou University of Science and Technology, 
on the basis of ipsr provided by Fei Hou and Chiyu Wang, Institute of Software, Chinese Academy of Sciences.

All rights reserved.

The codes can only be used for academic purpose, but cannot be used for commercial
purpose without written permission.

Redistribution and use in source for academic purpose, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include "kdtree.h"
#include "utility.h"
#include "PoissonRecon.h"
#include "Delaunay.h"
#include <set>
#include <map>
#include <Eigen/Dense>
#include <ctime>

#define PI 3.1415926
using namespace std;
typedef double REAL;
const unsigned int DIM = 3U;
Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
int iters = 10;
double pointweight = 10;
int depth = 10;

void save_points_normals(char* fname, vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	XForm<REAL, DIM + 1>& iXForm = XForm<REAL, DIM + 1>::Identity()) {
	FILE* fp = fopen(fname, "w");
	size_t n = points_normals.size();
	fprintf(fp, "ply\n");
	fprintf(fp, "format ascii 1.0\ncomment VCGLIB generated\n");
	fprintf(fp, "element vertex ");
	fprintf(fp, "%zd\n", n);
	fprintf(fp, "property float x\n");
	fprintf(fp, "property float y\n");
	fprintf(fp, "property float z\n");
	fprintf(fp, "property float nx\n");
	fprintf(fp, "property float ny\n");
	fprintf(fp, "property float nz\n");
	fprintf(fp, "element face 0\n");
	fprintf(fp, "property list uchar int vertex_indices\n");
	fprintf(fp, "end_header\n");
	for (int i = 0; i < n; i++) {
		Point<REAL, DIM> pt = iXForm * points_normals[i].first;
		Point<REAL, DIM> normal = points_normals[i].second.normal;
		fprintf(fp, "%f %f %f %f %f %f\n", pt[0], pt[1], pt[2], normal[0], normal[1], normal[2]);
	}
	fclose(fp);
}


void ransac_normal(vector<Point<REAL, DIM>>& points, size_t max_iters, double tol, Point<REAL, DIM>& ref_nrm,
	Point<REAL, DIM>& normal0, Point<REAL, DIM>& pt, size_t num0, double& min_dist) {

	size_t num = points.size();
	Point<REAL, DIM> center, normal, pt0, pt1, pt2, dir1, dir2;
	min_dist = 1e10;
	for (size_t iters = 0; iters < max_iters; iters++) {
		size_t id0 = rand() % num;
		pt0 = points[id0];
		center = pt0;
		size_t id1 = rand() % num;
		while (id0 == id1) 	id1 = rand() % num;
		pt1 = points[id1];
		dir1 = pt1 - pt0;
		center += pt1;
		size_t id2 = rand() % num;
		while (id2 == id1 || id2 == id0) id2 = rand() % num;
		pt2 = points[id2];
		dir2 = pt2 - pt0;
		center += pt2;
		center /= 3.0;
		Normalize(dir1), Normalize(dir2);
		normal = CrossProduct(dir1, dir2);
		Normalize(normal);

		size_t n = 3;
		for (size_t i = 0; i < num; i++) {
			if (i == id0 || i == id1 || i == id2) continue;
			pt1 = points[i];
			dir1 = pt1 - center;
			double dist = Point<REAL, DIM>::Dot(dir1, normal);
			if (fabs(dist) < tol) {
				n++;
			}
		}

		if (n >= num0) {
			Point<REAL, DIM> dev = pt - center;
			double dist = fabs(Point<REAL, DIM>::Dot(dev, normal));
			if (dist < min_dist) {
				normal0 = normal, min_dist = dist;
			}
		}
	}
	if (Point<REAL, DIM>::Dot(ref_nrm, normal0) < 0) normal0 *= -1;
}


bool reset_patches(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	vector<vector<int>>& sampleNbrs, vector<bool>& b_oriented, double tol0, double tol1, double tol2,
	vector<int>& orient_flag, vector< Point<REAL, DIM> >& ref_centers, vector<int>& curved_type, int t) {

	bool b_reset = false;
	size_t num = points_normals.size();
	vector< int > arr_oriented;

	size_t n = min(sampleNbrs[0].size(), t);
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>> patch1, patch2, patch3;
	for (int i = 0; i < num; i++) {
		if (orient_flag[i]) continue;  // 1邻域点要么全部已经定向，要么全部都没有定向
		queue<int> seeds;  // oriented
		if (b_oriented[i]) {
			seeds.push(i);
			while (!seeds.empty()) {
				int id0 = seeds.front(), flag = 1;
				Point<REAL, DIM> pt0 = points_normals[id0].first, nrm0 = points_normals[id0].second.normal;
				for (size_t j = 1; j < n; j++) {
					int id1 = sampleNbrs[id0][j];
					if (!b_oriented[id1]) {
						Point<REAL, DIM> pt1 = points_normals[id1].first, nrm1 = points_normals[id1].second.normal;
						bool b_updated = false;
						if (curved_type[id0] && j < 7) {
							Point<REAL, DIM> dev0 = ref_centers[id0] - pt0;  Normalize(dev0);
							double d0 = Point<REAL, DIM>::Dot(dev0, nrm0);
							Point<REAL, DIM> dev1 = ref_centers[id0] - pt1;  Normalize(dev1);
							double d4 = Point<REAL, DIM>::Dot(dev1, nrm1);
							if (fabs(d0) > tol0 && fabs(d4) > tol0) {
								if (d0 * d4 < 0) points_normals[id1].second.normal *= -1;
								b_updated = true;
								patch1.push_back(points_normals[id1]);
							}
						}
						else if (curved_type[id1] && j < 7) {
							Point<REAL, DIM> dev0 = ref_centers[id1] - pt0;  Normalize(dev0);
							double d0 = Point<REAL, DIM>::Dot(dev0, nrm0);
							Point<REAL, DIM> dev1 = ref_centers[id1] - pt1;  Normalize(dev1);
							double d4 = Point<REAL, DIM>::Dot(dev1, nrm1);
							if (fabs(d0) > tol0 && fabs(d4) > tol0) {
								if (d0 * d4 < 0) points_normals[id1].second.normal *= -1;
								b_updated = true;
								patch1.push_back(points_normals[id1]);
							}
						}
						else if (!b_updated) {
							Point<REAL, DIM> dev = pt1 - pt0;	Normalize(dev);
							double d1 = Point<REAL, DIM>::Dot(nrm0, nrm1), d2 = Point<REAL, DIM>::Dot(nrm0, dev), 
								d3 = Point<REAL, DIM>::Dot(nrm1, dev);
							if (fabs(d1) > tol2 && fabs(d2) < tol1 && fabs(d3) < tol1) {
								if (d1 < 0.0) points_normals[id1].second.normal *= -1;
								b_updated = true;
							}
						}
						if (b_updated) {
							b_oriented[id1] = true;
							b_reset = true;
							arr_oriented.push_back(id1);
							if (orient_flag[id1] != 1) seeds.push(id1);
						}
						else flag = 0;
					}
				}
				seeds.pop();
				orient_flag[id0] = flag;
			}
		}
	//}
		else {
			Point<REAL, DIM> pt0 = points_normals[i].first, nrm0 = points_normals[i].second.normal;
			int flag = -1;
			for (size_t j = 1; j < n; j++) {
				int id1 = sampleNbrs[i][j];
				if (b_oriented[id1]) {
					if (orient_flag[id1] != 1) seeds.push(id1);
					Point<REAL, DIM> pt1 = points_normals[id1].first, nrm1 = points_normals[id1].second.normal;
					bool b_updated = false;
					if (curved_type[id1] && j < 7) {
						Point<REAL, DIM> dev0 = ref_centers[id1] - pt0;   Normalize(dev0);
						double d0 = Point<REAL, DIM>::Dot(nrm0, dev0);
						Point<REAL, DIM> dev1 = ref_centers[id1] - pt1;   Normalize(dev1);
						double d4 = Point<REAL, DIM>::Dot(nrm1, dev1);
						if (fabs(d0) > tol0 && fabs(d4) > tol0) {
							if (d0 * d4 < 0) points_normals[i].second.normal *= -1;
							b_updated = true;
							patch2.push_back(points_normals[i]);
						}
					}
					else if (curved_type[i] && j < 7) {
						Point<REAL, DIM> dev0 = ref_centers[i] - pt0;   Normalize(dev0);
						double d0 = Point<REAL, DIM>::Dot(nrm0, dev0);
						Point<REAL, DIM> dev1 = ref_centers[i] - pt1;   Normalize(dev1);
						double d4 = Point<REAL, DIM>::Dot(nrm1, dev1);
						if (fabs(d0) > tol0 && fabs(d4) > tol0) {
							if (d0 * d4 < 0) points_normals[i].second.normal *= -1;
							b_updated = true;
							patch2.push_back(points_normals[i]);
						}
					}
					else if (!b_updated) {
						Point<REAL, DIM> dev = pt1 - pt0;	Normalize(dev);
						double d1 = Point<REAL, DIM>::Dot(nrm0, nrm1), d2 = Point<REAL, DIM>::Dot(nrm0, dev), d3 = Point<REAL, DIM>::Dot(nrm1, dev);
						if (fabs(d1) > tol2 && fabs(d2) < tol1 && fabs(d3) < tol1) {
							if (d1 < 0.0) points_normals[i].second.normal *= -1;
							b_updated = true;
						}
					}
					if (b_updated) {
						flag = 0;
						b_oriented[i] = true;
						b_reset = true;
						arr_oriented.push_back(i);
						seeds.push(i);
						break;
					}
				}
			}

			orient_flag[i] = flag;
			while (!seeds.empty()) {
				int id0 = seeds.front(), flag = 1;
				Point<REAL, DIM> pt0 = points_normals[id0].first, nrm0 = points_normals[id0].second.normal;
				for (size_t j = 1; j < n; j++) {
					int id1 = sampleNbrs[id0][j];
					if (!b_oriented[id1]) {
						Point<REAL, DIM> pt1 = points_normals[id1].first, nrm1 = points_normals[id1].second.normal;
						bool b_updated = false;
						if (curved_type[id0] && j < 7) {
							Point<REAL, DIM> dev0 = ref_centers[id0] - pt0;  Normalize(dev0);
							double d0 = Point<REAL, DIM>::Dot(dev0, nrm0);
							Point<REAL, DIM> dev1 = ref_centers[id0] - pt1;  Normalize(dev1);
							double d4 = Point<REAL, DIM>::Dot(dev1, nrm1);
							if (fabs(d0) > tol0 && fabs(d4) > tol0) {
								if (d0 * d4 < 0) points_normals[id1].second.normal *= -1;
								b_updated = true;
								patch3.push_back(points_normals[id1]);
							}
						}
						else if (curved_type[id1] && j < 7) {
							Point<REAL, DIM> dev0 = ref_centers[id1] - pt0;  Normalize(dev0);
							double d0 = Point<REAL, DIM>::Dot(dev0, nrm0);
							Point<REAL, DIM> dev1 = ref_centers[id1] - pt1;  Normalize(dev1);
							double d4 = Point<REAL, DIM>::Dot(dev1, nrm1);
							if (fabs(d0) > tol0 && fabs(d4) > tol0) {
								if (d0 * d4 < 0) points_normals[id1].second.normal *= -1;
								b_updated = true;
								patch3.push_back(points_normals[id1]);
							}
						}
						else if (!b_updated) {
							Point<REAL, DIM> dev = pt1 - pt0;	Normalize(dev);
							double d1 = Point<REAL, DIM>::Dot(nrm0, nrm1), d2 = Point<REAL, DIM>::Dot(nrm0, dev),
								d3 = Point<REAL, DIM>::Dot(nrm1, dev);
							if (fabs(d1) > tol2 && fabs(d2) < tol1 && fabs(d3) < tol1) {
								if (d1 < 0.0) points_normals[id1].second.normal *= -1;
								b_updated = true;
							}
						}
						if (b_updated) {
							b_oriented[id1] = true;
							b_reset = true;
							arr_oriented.push_back(id1);
						}
						else flag = 0;
					}
				}
				seeds.pop();
				orient_flag[id0] = flag;
			}
		}
	}

	return b_reset;
}


Point<REAL, DIM> normal_pca(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals, vector < int >& sub_structure,
	Point<REAL, DIM>& centroid, size_t num) {

	centroid = Point<REAL, DIM>(0, 0, 0);
	Eigen::MatrixXd matT(num, 3);
	for (size_t i = 0; i < num; ++i) {
		int id = sub_structure[i];
		centroid += points_normals[id].first;
		matT(i, 0) = points_normals[id].first[0];
		matT(i, 1) = points_normals[id].first[1];
		matT(i, 2) = points_normals[id].first[2];
	}
	centroid /= num;

	for (size_t i = 0; i < num; i++)
	{
		matT(i, 0) -= centroid[0];
		matT(i, 1) -= centroid[1];
		matT(i, 2) -= centroid[2];
	}

	Eigen::MatrixXd matA = matT.transpose() * matT;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matA);
	Eigen::MatrixXd eigenVectors = eigensolver.eigenvectors();
	return Point<REAL, DIM> (eigenVectors(0, 0), eigenVectors(1, 0), eigenVectors(2, 0));
}


// ransac
Point<REAL, DIM> normal_ransac(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	vector < int >& sub_structure, size_t n, double tol) {

	size_t n0 = 0, num = min(sub_structure.size(), n);
	Point<REAL, DIM> center, center0, normal, normal0, pt0, pt1, pt2, dir1, dir2;
	for (size_t iters = 0; iters < 10; iters++) {
		size_t id0 = rand() % num;
		pt0 = points_normals[sub_structure[id0]].first;
		center = pt0;
		size_t id1 = rand() % num;
		while (id0 == id1) 	id1 = rand() % num;
		pt1 = points_normals[sub_structure[id1]].first;
		dir1 = pt1 - pt0;
		center += pt1;
		size_t id2 = rand() % num;
		while (id2 == id1 || id2 == id0) id2 = rand() % num;
		pt2 = points_normals[sub_structure[id2]].first;
		dir2 = pt2 - pt0;
		center += pt2;
		center /= 3.0;
		Normalize(dir1), Normalize(dir2);
		normal = CrossProduct(dir1, dir2);
		Normalize(normal);

		size_t n = 3;
		for (size_t i = 0; i < num; i++) {
			if (i == id0 || i == id1 || i == id2) continue;
			pt1 = points_normals[sub_structure[i]].first;
			dir1 = pt1 - center;
			double dist = Point<REAL, DIM>::Dot(dir1, normal);
			if (fabs(dist) < tol) {
				n++;
			}
			else if (i == 0) {
				break;
			}
		}

		if (n > n0) {
			n0 = n, normal0 = normal, center0 = center;
		}
		if (n0 == num) break;
	}
	return normal0;
}

void normals_smooth(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	vector<vector<int>>& sampleNbrs, size_t n, double tol) {

	size_t num = min(sampleNbrs[0].size(), n);
	for (size_t i = 0; i < points_normals.size(); i++) {
		vector< Point<REAL, DIM> > arr_nrmls;
		Point<REAL, DIM> nrm0 = points_normals[i].second.normal;
		arr_nrmls.push_back(nrm0);
		for (size_t j = 1; j < num; j++) {
			Point<REAL, DIM> nrm1 = points_normals[sampleNbrs[i][j]].second.normal;
			bool is_merged = false;
			for (size_t k = 0; k < arr_nrmls.size(); k++) {
				Point<REAL, DIM> nrm2 = arr_nrmls[k];
				Normalize(nrm2);
				double cosin = Point<REAL, DIM>::Dot(nrm1, nrm2);
				if (fabs(cosin) > tol) {
					if (cosin > 0)	arr_nrmls[k] += nrm1;
					else			arr_nrmls[k] -= nrm1;
					is_merged = true;
					break;
				}
			}
			if (!is_merged)	arr_nrmls.push_back(nrm1);
		}
		double max_len = 0.0;
		size_t idx = 0;
		for (size_t k = 0; k < arr_nrmls.size(); k++) {
			double len = Length(arr_nrmls[k]);
			if (len > max_len) {
				max_len = len;
				idx = k;
			}
		}
		Point<REAL, DIM> nrm = arr_nrmls[idx];
		Normalize(nrm);
		if (fabs(Point<REAL, DIM>::Dot(nrm, nrm0)) < tol) {
			points_normals[i].second.normal = nrm;
		}
	}
}


void calc_normals(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals, vector<vector<int>>& sampleNbrs, 
	vector<double>& arr_lens, size_t num, double ave_len, double f1=0.3, double f2=1.5) {

	Point<REAL, DIM> centroid(0, 0, 0);
	size_t n = points_normals.size(), s = min(num, sampleNbrs[0].size());
	double tol1 = ave_len * f1, tol2 = ave_len * f2;
	for (size_t i = 0; i < n; i++) {
		if (arr_lens[i] < tol2)
			points_normals[i].second.normal = normal_pca(points_normals, sampleNbrs[i], centroid, s);
		else 
			points_normals[i].second.normal = normal_ransac(points_normals, sampleNbrs[i], s, tol1);
	}
	normals_smooth(points_normals, sampleNbrs, s, 0.6);
}


bool is_visible (Point<REAL, DIM>& ref_pnt, double rd, double& md,
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals) {

	md = 1e10;
	size_t s = 10;
	for (size_t i = 0; i < points_normals.size()/s; i++) {
		int id = i * s + rand() % s;
		Point<REAL, DIM> pt = points_normals[id].first;
		double d = Distance(ref_pnt, pt);
		if (d < md) md = d;
		if (d < rd) return false;
	}
	return true;
}



bool cacl_center(vector < Point<REAL, DIM>>& patches, Point<REAL, DIM>& center, double tol1, double tol2) {

	size_t n = patches.size();
	Eigen::MatrixXd matA(3, 3), matB(3, 1), matT(3, 1);
	matA = Eigen::MatrixXd::Zero(3, 3);
	center = Point<REAL, DIM>(0, 0, 0);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			Point<REAL, DIM> dev = patches[i] - patches[j];
			matT(0, 0) = dev[0], matT(1, 0) = dev[1], matT(2, 0) = dev[2];
			matA += matT * matT.transpose();
			double d = Point<REAL, DIM>::Dot(patches[i], patches[i]) -
				Point<REAL, DIM>::Dot(patches[j], patches[j]);
			center += dev * d;
		}
	}
	center *= 0.5;
	matB(0, 0) = center[0], matB(1, 0) = center[1], matB(2, 0) = center[2];
	matB = matA.inverse() * matB;
	center[0] = matB(0, 0), center[1] = matB(1, 0), center[2] = matB(2, 0);

	double r = 0.0;
	vector<double> arr_radius(n);
	for (size_t i = 0; i < n; i++) {
		double d = Distance(center, patches[i]);
		arr_radius[i] = d;
		r += d;
	}
	r /= n;
	if (r > tol2) {
		return false;
	}

	double std_var = 0.0;
	for (size_t i = 0; i < n; i++) {
		double d = arr_radius[i] - r;
		std_var += d * d;
	}
	std_var = sqrt(std_var / (n - 1));
	if (std_var / r < tol1)	return true;
	return false;
}



void orient_normals(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	vector<vector<int>>& sampleNbrs, double diameter, double r,
	double ave_len, size_t n2, XForm<REAL, DIM + 1>& iXForm) {
	// n1-neighbor used to calc normals; n2-neighbor used to propagate normal orientation

	size_t n = points_normals.size();
	double rd = diameter * 0.998;
	vector<bool> b_oriented(n, false);  
	// 先对部分法矢作初步定向，作为后续法矢定向的guidance
	vector< Point<REAL, DIM> > ref_centers(n, Point<REAL, DIM>(0, 0, 0));
	vector<int> curved_type(n, 0);
	Point < REAL, DIM> center, normal;
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>> patch0;


//#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		Point<REAL, DIM> normal = points_normals[i].second.normal;
		Point<REAL, DIM> pt = points_normals[i].first;
		Point<REAL, DIM> ref1 = pt - normal * diameter, ref2 = pt + normal * diameter;

		// 利用距离推测可见性从而定向法矢
		double d1, d2;
		bool b1 = is_visible(ref1, rd, d1, points_normals), b2 = is_visible(ref2, rd, d2, points_normals);
		if (b1 || b2) {
			b_oriented[i] = true;
			if (b1 && (b2 && d1 > d2) || !b2) {
				normal *= -1;
				points_normals[i].second.normal = normal;
			}
		}
		if (b_oriented[i])	patch0.push_back(points_normals[i]);
	}

	patch0.clear();
	size_t m = sampleNbrs[0].size();
	vector<Point<REAL, DIM>> patches(m, Point<REAL, DIM>(0, 0, 0));
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>> patches11;
	double len1 = ave_len * 1.4, len2 = ave_len *1.5;
	for (int i = 0; i < n; i++) {
		for (size_t j = 0; j < m; j++) {
			patches[j] = points_normals[sampleNbrs[i][j]].first;
		}
		if (cacl_center(patches, center, r, len1)) {
			curved_type[i] = 1;
			pair<Point<REAL, DIM>, Normal<REAL, DIM>> pn = points_normals[i];
			normal = pn.first - center;
			Normalize(normal);
			ref_centers[i] = pn.first - normal * len2;
			pn.second = -normal;
			patch0.push_back(pn);
		}
	}
	// 以已定向顶点为“种子”，用区域增长法设置其他顶点的法矢
	double t0 = 0.9, t1 = 0.12, t2 = 0.98;
	vector<int> orient_flag(n, 0);  // 1：邻域完全定向；-1：邻域完全未定向；0；部分定向
	size_t iter = 0, max_iter = 22;
	while (reset_patches(points_normals, sampleNbrs, b_oriented, t0, t1, t2,
		orient_flag, ref_centers, curved_type, n2)) {

		iter++;
		t0 /= 1.02;  
		t1 *= 1.07;
		t2 /= 1.07;
		patch0.clear();
		for (int i = 0; i < n; i++) {
			if (b_oriented[i]) patch0.push_back(points_normals[i]);
		}

		printf("%zd:%zd / %zd\n", iter, patch0.size(), n);
		if (iter == max_iter) break;
	}

	for (int i = 0; i < n; i++) {
		if (b_oriented[i]) continue;
		Point<REAL, DIM> normal0 = points_normals[i].second.normal;
		Point<REAL, DIM> normal1(0, 0, 0);
		for (int j = 1; j < 7; j++) {
			normal1 += points_normals[sampleNbrs[i][j]].second.normal;
		}
		normal1 /= 6;
		double d = Point<REAL, DIM>::Dot(normal0, normal1);
		if (fabs(d) > 0.25) {
			if (d < 0)	points_normals[i].second.normal = -normal0;
		}
	}
}



#include<cstdlib>
bool is_noisy(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	vector<vector<int>> sampleNbrs, double tol) {

	double d = 0.0;
	size_t num = points_normals.size(), s = sampleNbrs[0].size() - 1;
	for (size_t i = 0; i < 100; i++) {
		size_t r = rand() % num;
		Point<REAL, DIM> pt0 = points_normals[r].first, nrml = points_normals[r].second.normal;
		for (size_t j = 1; j <= s; j++) {
			Point<REAL, DIM> pt1 = points_normals[sampleNbrs[r][j]].first;
			Point<REAL, DIM> nrml1 = points_normals[sampleNbrs[r][j]].second.normal;
			Point<REAL, DIM> dev = pt1 - pt0;
			Normalize(dev);
			d += fabs(Point<REAL, DIM>::Dot(dev, nrml)) + fabs(Point<REAL, DIM>::Dot(dev, nrml1));
		}
	}
	d /= (200 * s);
	bool b = d > tol ? true : false;
	return b;
}

void filter_normals(vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>>& points_normals,
	vector<vector<int>>& sampleNbrs, double sigma, size_t nbr_size, size_t iter_num) {

	size_t sample_num = sampleNbrs.size();

	// 设置距离权
	vector<vector<double>> distance_wts;
	distance_wts.resize(sample_num);
	vector < Point<REAL, DIM>> temp_nrmls;
	temp_nrmls.resize(sample_num);
	for (size_t i = 0; i < sample_num; i++) {
		distance_wts[i].reserve(nbr_size);
		Point<REAL, DIM> pt0 = points_normals[i].first;
		for (size_t j = 1; j <= nbr_size; j++) {
			Point<REAL, DIM> pt1 = points_normals[sampleNbrs[i][j]].first;
			Point<REAL, DIM> dev = pt0 - pt1;
			double d = Point<REAL, DIM>::Dot(dev, dev);
			distance_wts[i].push_back(exp(-d / sigma));
		}
	}

	for (size_t iter = 0; iter < iter_num; iter++) {
		for (size_t i = 0; i < sample_num; i++) {
			Point<REAL, DIM> nrml = points_normals[i].second.normal * 0.4;
			for (size_t j = 1; j <= nbr_size; j++) {
				int k = sampleNbrs[i][j];
				double w = distance_wts[i][j - 1];
				double nrml_wt = Point<REAL, DIM>::Dot(points_normals[i].second.normal, points_normals[k].second.normal);
				if (nrml_wt < 0) w *= -exp((nrml_wt * nrml_wt - 1) * 9);
				else w *= exp((nrml_wt * nrml_wt - 1) * 9);
				nrml += points_normals[k].second.normal * w;
			}
			Normalize(nrml);
			temp_nrmls[i] = nrml;
		}
		for (size_t i = 0; i < sample_num; i++) {
			points_normals[i].second.normal = temp_nrmls[i];
		}
	}
}



int main(int argc, char* argv[])
{
	string input_name = "examples/angel.ply", output_name1 = "examples/angel_gnop.ply";
	char output_name2[] = "examples/angel_gnop_o.ply";
	size_t k1 = 10, k2 = 7;          // important parameter  // 噪声大时应增大k1
	double t2 = 0.14;

	clock_t start, end;
	start = clock();
	vector<pair<Point<REAL, DIM>, Normal<REAL, DIM>>> points_normals;
	ply_reader<REAL, DIM>(input_name, points_normals);

	string command = "PoissonRecon --in i.ply --out o.ply --bType 2 --depth " + to_string(depth) + " --pointWeight " + to_string(pointweight);
	vector<string> cmd = split(command);
	vector<char*> argv_str(cmd.size());
	for (size_t i = 0; i < cmd.size(); ++i)
		argv_str[i] = &cmd[i][0];

	XForm<REAL, DIM + 1> iXForm;
	vector<double> weight_samples;
	double ave_len = 0.0, xm = 1e10, ym = 1e10, zm = 1e10, xM = -1e10, yM = -1e10, zM = -1e10, diameter;
	points_normals = SamplePoints<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, iXForm, &weight_samples);
	size_t samples_num = points_normals.size();
	vector<vector<int>> sampleNbrs;
	sampleNbrs.resize(samples_num);
	kdt::KDTree<kdt::KDTreePoint> tree;
	vector<kdt::KDTreePoint> vertices;
	vertices.reserve(samples_num);
	Point<REAL, DIM> center(0, 0, 0);
	for (size_t i = 0; i < samples_num; ++i) {
		center += points_normals[i].first;
		array<double, 3> p{ points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2] };
		if (xm > points_normals[i].first[0])		xm = points_normals[i].first[0];
		else if (xM < points_normals[i].first[0])	xM = points_normals[i].first[0];
		if (ym > points_normals[i].first[1])		ym = points_normals[i].first[1];
		else if (yM < points_normals[i].first[1])   yM = points_normals[i].first[1];
		if (zm > points_normals[i].first[2])		zm = points_normals[i].first[2];
		else if (zM < points_normals[i].first[2])	zM = points_normals[i].first[2];
		vertices.push_back(kdt::KDTreePoint(p));
	}
	tree.build(vertices);
	center /= samples_num;
	Point<REAL, DIM> corner_bl = Point<REAL, DIM>(xm, ym, zm), corner_ur = Point<REAL, DIM>(xM, yM, zM);

	diameter = Distance(corner_bl, corner_ur);
	pair<vector<Point<REAL, DIM>>, vector<vector<int>>> mesh;
	vector<double> arr_lens(samples_num, 0);
	for (size_t i = 0; i < samples_num; ++i) {
		sampleNbrs[i] = tree.knnSearch(vertices[i], k1);
		size_t id0 = sampleNbrs[i][0];
		size_t id1 = sampleNbrs[i][6];
		arr_lens[i] = Distance(points_normals[id0].first, points_normals[id1].first);
	}
	size_t num = static_cast<size_t>(samples_num / 50);
	for (size_t i = 0; i < num; i++) {
		ave_len += arr_lens[i*50];
	}
	ave_len /= num;

	calc_normals(points_normals, sampleNbrs, arr_lens, k1, ave_len*5);  // 非均匀时应减小p5,参考值：ave_len*5
//	save_points_normals("examples/angel_unoriented.ply", points_normals, iXForm);
	if (is_noisy(points_normals, sampleNbrs, 0.15)) {
		filter_normals(points_normals, sampleNbrs, ave_len * ave_len, 6, 5);
		printf_s("filter_normals\n");
	}
	orient_normals(points_normals, sampleNbrs, diameter, t2, ave_len, k2, iXForm); // 噪声大时应增大p4
	end = clock();
	double sec = (end - start) / 1000.0;
	printf_s("Running time %f s\n", sec);
	mesh = poisson_reconstruction<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, &weight_samples);
	end = clock();
	sec = (end - start) / 1000.0;
	printf_s("Running time %f s\n", sec);

	end = clock();
	sec = (end - start) / 1000.0;
	printf_s("Running time %f s\n", sec);
	save_points_normals(output_name2 , points_normals, iXForm);
	output_ply(output_name1, mesh, iXForm);
	return 0;
}
