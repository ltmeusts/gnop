#pragma once
#include<math.h>
#include<vector>
#include<algorithm>

#include"../MeshLib/core/Mesh/Vertex.h"
#include"../MeshLib/core/Mesh/Edge.h"
#include"../MeshLib/core/Mesh/Face.h"
#include"../MeshLib/core/Mesh/HalfEdge.h"
#include"../MeshLib/core/Mesh/BaseMesh.h"
#include"../MeshLib/core/Mesh/iterators.h"


using namespace MeshLib;
using namespace std;

typedef CBaseMesh<CVertex, CEdge, CFace, CHalfEdge> Mesh;
typedef VertexEdgeIterator<CVertex, CEdge, CFace, CHalfEdge> veIter;
typedef VertexFaceIterator<CVertex, CEdge, CFace, CHalfEdge> vfIter;
typedef FaceVertexIterator<CVertex, CEdge, CFace, CHalfEdge> fvIter;
typedef FaceEdgeIterator<CVertex, CEdge, CFace, CHalfEdge> feIter;
typedef FaceHalfedgeIterator<CVertex, CEdge, CFace, CHalfEdge> fhIter;
typedef MeshVertexIterator<CVertex, CEdge, CFace, CHalfEdge> mvIter;
typedef MeshEdgeIterator<CVertex, CEdge, CFace, CHalfEdge> meIter;
typedef MeshFaceIterator<CVertex, CEdge, CFace, CHalfEdge> mfIter;


class MyPoint2 {
public:
	double x;
	double y;
	MyPoint2() {}
	MyPoint2(double x, double y) {
		this->x = x;
		this->y = y;
	}
};

int delaunayVertexId = 0;  ///????????
int delaunayFaceId = 0;

//�������򣬴�С������x�����С������x������ͬ����y�����С����
void insertSort(vector<CVertex*>& allVertex) {
	int i, j;
	CVertex* tmp = NULL;
	for (i = 0; i < allVertex.size(); i++) {
		tmp = allVertex[i];
		j = i - 1;
		while (j >= 0 && allVertex[j]->point().coord()[0] >= tmp->point().coord()[0]) {
			if (allVertex[j]->point().coord()[0] == tmp->point().coord()[0]
				&& allVertex[j]->point().coord()[1] <= tmp->point().coord()[1]) {
				break;
			}
			allVertex[j + 1] = allVertex[j];
			j--;
		}
		allVertex[j + 1] = tmp;
	}
}

//������������0��ʾ��oa��obΪ��ʱ����ת
double cross(CVertex* o, CVertex* a, CVertex* b) {
	return (a->point().coord()[0] - o->point().coord()[0])
		* (b->point().coord()[1] - o->point().coord()[1])
		- (a->point().coord()[1] - o->point().coord()[1])
		* (b->point().coord()[0] - o->point().coord()[0]);
}

//Andrew's Monotone Chain
void ConvexHull(vector<CVertex*>& allVertex, vector<CVertex*>& vecConvexHull) {
	//��Vertex�Ķ���id����
	insertSort(allVertex);

	//��¼��Χ��ĸ���
	int m = 0;
	size_t size = allVertex.size();

	//�����ң��ó��°벿�ְ�Χ��
	for (int i = 0; i < size; i++) {
		while (m >= 2 && cross(vecConvexHull[m - 2], vecConvexHull[m - 1], allVertex[i]) < 0) {
			vecConvexHull.pop_back();
			m--;
		}
		vecConvexHull.push_back(allVertex[i]);
		m++;
	}

	//���ҵ��󣬵ó��ϰ벿�ְ�Χ�㣬��洢��һ�����id
	for (size_t i = size - 2; i >= 0; i--) {
		while (cross(vecConvexHull[m - 2], vecConvexHull[m - 1], allVertex[i]) < 0) {
			vecConvexHull.pop_back();
			m--;
		}
		vecConvexHull.push_back(allVertex[i]);
		m++;
	}

	vecConvexHull.pop_back();
}

//����Բ�ģ�����Ϊ������
MyPoint2 getCenter(vector<CVertex*> vec) {
	//������Բ��Բ��
	double x0 = vec[0]->point().coord()[0];
	double y0 = vec[0]->point().coord()[1];
	double x1 = vec[1]->point().coord()[0];
	double y1 = vec[1]->point().coord()[1];
	double x2 = vec[2]->point().coord()[0];
	double y2 = vec[2]->point().coord()[1];

	double a0 = (x0 + x1) / 2;
	double b0 = (y0 + y1) / 2;
	double a1 = (x2 + x1) / 2;
	double b1 = (y2 + y1) / 2;

	//��������
	double x;
	double y;


	if (y0 != y1 && y1 != y2) {
		double k1 = -(x0 - x1) / (y0 - y1);
		double c1 = b0 - k1 * a0;
		double k2 = -(x1 - x2) / (y1 - y2);
		double c2 = b1 - k2 * a1;

		x = (c2 - c1) / (k1 - k2);
		y = k1 * x + c1;
	}
	else if (y0 == y1) {
		x = a0;
		double k2 = -(x1 - x2) / (y1 - y2);
		double c2 = b1 - k2 * a1;
		y = k2 * x + c2;
	}
	else {
		x = a1;
		double k1 = -(x0 - x1) / (y0 - y1);
		double c1 = b0 - k1 * a0;
		y = k1 * x + c1;
	}

	MyPoint2 point2(x, y);

	return point2;
}

//����Բ�ģ�����Ϊһ����
MyPoint2 getTriCenter(CFace* pFace) {
	vector<CVertex*> vec;
	for (fvIter iter(pFace); !iter.end(); iter++) {
		vec.push_back(*iter);
	}
	return getCenter(vec);
}

//��Face����Ӧ��Բ��Vertex����mpFaceCenter��
void addMpFaceCenter(map<CFace*, MyPoint2>& mpFaceCenter, CFace* pFace) {
	MyPoint2 point2 = getTriCenter(pFace);
	mpFaceCenter.insert(pair<CFace*, MyPoint2>(pFace, point2));
}

//��Face����Ӧ��Բ��Vertex��mpFaceCenter��ɾ��
void delMpFaceCenter(map<CFace*, MyPoint2>& mpFaceCenter, CFace* pFace) {
	mpFaceCenter.erase(pFace);
}

//Euclidean Distance��ƽ��
double eDistance(double x1, double y1, double x2, double y2) {
	return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

//�ж�һ�����Ƿ���һ������������Բ��
bool isCircleContain(map<CFace*, MyPoint2>& mpFaceCenter, CFace* pFace, CVertex* pVertex) {
	MyPoint2 point2 = mpFaceCenter.at(pFace);
	fvIter iter(pFace);

	double r = eDistance(point2.x, point2.y, (*iter)->point().coord()[0], (*iter)->point().coord()[1]);
	double distance = eDistance(point2.x, point2.y, pVertex->point().coord()[0], pVertex->point().coord()[1]);

	return r > distance;
}

//����ConvexHull���ɵĵ㼯��ʼ��Delaunayͼ
void initDelaunayMesh(const vector<CVertex*>& allVertex, const vector<CVertex*>& vecConvexHull, 
	Mesh* pDelaunayMesh, map<CFace*, MyPoint2>& mpFaceCenter) {
	vector<CVertex*> vecTmp = vecConvexHull;
	bool isContainOther = false;
	while (vecTmp.size() > 3) {
		size_t size = vecTmp.size();
		for (int i = 0; i < size; i++) {
			isContainOther = false;
			vector<CVertex*> fv;
			fv.push_back(vecTmp[i]);
			fv.push_back(vecTmp[(i + 1) % size]);
			fv.push_back(vecTmp[(i + 2) % size]);

			MyPoint2 center = getCenter(fv);
			double r = eDistance(center.x, center.y, fv[0]->point().coord()[0], fv[0]->point().coord()[1]);

			for (int j = 0; j < vecConvexHull.size(); j++) {
				CVertex* pVertexTmp = vecConvexHull[j];
				if (pVertexTmp != fv[0] && pVertexTmp != fv[1] && pVertexTmp != fv[2]
					&& r > eDistance(center.x, center.y, pVertexTmp->point().coord()[0], pVertexTmp->point().coord()[1])) {
					isContainOther = true;
					break;
				}
			}
			if (!isContainOther) {
				mpFaceCenter.insert(pair<CFace*, MyPoint2>((*pDelaunayMesh).createFace(fv, delaunayFaceId++), center));
				vecTmp.erase(vecTmp.begin() + (i + 1) % size);
				break;
			}
		}
	}

	addMpFaceCenter(mpFaceCenter, (*pDelaunayMesh).createFace(vecTmp, delaunayFaceId++));
}

//�ж��Ƿ�ΪDelaunayTriangle
bool isDelaunayTriangle(Mesh* pDelaunayMesh, map<CFace*, MyPoint2>& mpFaceCenter, vector<CVertex*>& allVertex) {
	for (mfIter iter(pDelaunayMesh); !iter.end(); iter++) {
		fvIter iter2(*iter);
		CVertex* v1 = *iter2;
		iter2++;
		CVertex* v2 = *iter2;
		iter2++;
		CVertex* v3 = *iter2;
		for (int i = 0; i < allVertex.size(); i++) {
			if (v1 != allVertex[i] && v2 != allVertex[i] && v3 != allVertex[i] && 
				isCircleContain(mpFaceCenter, *iter, allVertex[i])) {
				return false;
			}
		}
	}
	return true;
}

//Delaunay Triangulation
void DelaunayTriangulation(vector<CVertex*>& allVertex, Mesh* pDelaunayMesh) {

	// ��ȡ͹�߽�
	vector<CVertex*> vecConvexHull;
	ConvexHull(allVertex, vecConvexHull);

	//��ʼ��toBeAdd
	vector<CVertex*> toBeAdd;
	for (int i = 0; i < allVertex.size(); i++) {
		size_t num = count(vecConvexHull.begin(), vecConvexHull.end(), allVertex[i]);
		if (num == 0) {
			toBeAdd.push_back(allVertex[i]);
		}
	}

	//�ð�Χ��������һ�㴴��������
	map<CFace*, MyPoint2> mpFaceCenter;
	initDelaunayMesh(allVertex, vecConvexHull, pDelaunayMesh, mpFaceCenter);  ////????????????????????

	//�����DelaunayMesh����ӵ�
	for (int n = 0; n < toBeAdd.size(); n++) {
		//ȡ���������Բ������ǰ���Face
		vector<CFace*> vecFace;
		for (mfIter iter(pDelaunayMesh); !iter.end(); iter++) {
			if (isCircleContain(mpFaceCenter, *iter, toBeAdd[n])) {
				vecFace.push_back(*iter);
			}
		}

		//ȡ����һ���õ�������Face�໥֮�䲻����ıߵĶ˵�,������ʱ��˳�����vecVerTmp��
		vector<CVertex*> vecVerTmp;
		for (int i = 0; i < vecFace.size(); i++) {
			CFace* f0 = vecFace[i];
			for (feIter iter(f0); !iter.end(); iter++) {
				CEdge* tmpEdge = *iter;
				CHalfEdge* h0 = tmpEdge->halfedge(0);
				CHalfEdge* h1 = tmpEdge->halfedge(1);
				assert(h0 != NULL);
				if (h1 == NULL) {
					vecVerTmp.push_back(h0->source());
					vecVerTmp.push_back(h0->target());
				}
				else {
					CFace* f1 = (f0 == (*pDelaunayMesh).edgeFace1(tmpEdge))
						? (*pDelaunayMesh).edgeFace2(tmpEdge)
						: (*pDelaunayMesh).edgeFace1(tmpEdge);
					vector<CFace*>::iterator tmpIter = find(vecFace.begin(), vecFace.end(), f1);
					if (tmpIter == vecFace.end()) {

						CHalfEdge* tmpHe = f0->halfedge();
						while (tmpHe != h0 && tmpHe != h1) {
							tmpHe = tmpHe->he_next();
						}
						vecVerTmp.push_back(tmpHe->source());
						vecVerTmp.push_back(tmpHe->target());
					}
				}
			}
		}

		//ɾ�����е�һ���õ���Face
		for (int i = 0; i < vecFace.size(); i++) {
			delMpFaceCenter(mpFaceCenter, vecFace[i]);
			(*pDelaunayMesh).deleteFace(vecFace[i]);
		}

		//ÿ��������ߵ������˵��뵱ǰ���γ��µ�Face
		for (int i = 0; i < vecVerTmp.size(); i += 2) {
			vector<CVertex*> tmp;
			tmp.push_back(toBeAdd[n]);
			tmp.push_back(vecVerTmp[i]);
			tmp.push_back(vecVerTmp[i + 1]);
			addMpFaceCenter(mpFaceCenter, (*pDelaunayMesh).createFace(tmp, delaunayFaceId++));
		}
	}
	pDelaunayMesh->labelBoundary();
}