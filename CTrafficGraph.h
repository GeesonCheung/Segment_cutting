#ifndef CTRAFICGRAPH_H_
#define CTRAFICGRAPH_H_

#include <map>
#include <vector>
using namespace std;

/*
    结构体：图模型中的点数据结构
*/
typedef struct stfPoint {
	double x;
	double y;
	unsigned int pid;
	stfPoint()
	{
		
	}
	stfPoint(float _x,float _y)
	{
		x = _x;
		y = _y;
	}
	stfPoint(int _id, float _x, float _y)
	{
		pid = _id;
		x = _x;
		y = _y;
	}
	bool operator < (const stfPoint &other) const
	{
		if (doubleEqual(this->x, other.x) && doubleEqual(this->y, other.y))
			return false;
		else {
			
			if (this->pid < other.pid)
				return true;
			else
				return false;
		}
	}

	bool operator == (const stfPoint &other) const
	{
		if (doubleEqual(this->x, other.x) && doubleEqual(this->y, other.y))
			return true;
		else
			return false;
	}

}stfPoint_t;

/*
    结构体：图模型中的边数据结构
*/
typedef struct stfEdge {
	stfPoint StartPoint;
	stfPoint EndPoint;
	float Weight;
	unsigned int eid;
	stfEdge()
	{
		
	}
	stfEdge(stfPoint_t p1,stfPoint_t p2)
	{
		StartPoint = p1;
		EndPoint = p2;
	}
	stfEdge(stfPoint_t p1, stfPoint_t p2,float weight)
	{
		StartPoint = p1;
		EndPoint = p2;
		Weight = weight;
	}
	
	bool operator < (const stfEdge& other) const
	{
		
		if ((StartPoint == other.StartPoint && EndPoint == other.EndPoint) ||
			(StartPoint == other.EndPoint && EndPoint == other.StartPoint))
			return  false;
		else {
			if (this->eid < other.eid)
				return true;
			else
				return false;
		}
	}

	bool operator == (const stfEdge &other) const
	{
		if ((StartPoint == other.StartPoint && EndPoint == other.EndPoint) ||
			(StartPoint == other.EndPoint && EndPoint == other.StartPoint))
			return true;
		else
			return false;
	}

}stfEdge_t;

typedef unsigned long uint64;

/*
  建图算法类
*/
class DLL_CLASS_API CTrafficGraph
{
public:
	CTrafficGraph(vector<stfEdge>& vecEdges);
    ~CTrafficGraph();

public:
	void DumpEdge();

public:
	/*
	   *Summary: 向已有线段集合中添加线段
		*Parameters:
		      @input: 
			 startPoint: 线段起点 
		         endPoint:   线段终点
		*Return: void.
	*/
	void InputEdge(stfPoint_t& startPoint, stfPoint_t& endPoint);

	/*
	   图模型中的点集合
	*/
	map < stfPoint_t, unsigned int >	fmapPoints;
	/*
	   图模型中的边集合
	*/
	map < stfEdge_t, unsigned int >		fmapEdges;

public:
	bool lineIntersection(stfPoint_t& startPoint, stfPoint_t& endPoint, stfEdge_t& edge, stfPoint_t& crossPoint);
	bool PointEqual(const stfPoint_t& Point1,const stfPoint_t& Point2);
	bool hasPoint(const stfPoint_t& Point, unsigned int& ID);
	bool hasEdge(const stfPoint_t& startPoint,const stfPoint_t& endPoint, unsigned int& ID);
	unsigned int  AddPoint( stfPoint_t& Point);
	unsigned int  AddEdge( stfPoint_t& startPoint, stfPoint_t& endPoint);
	double PointToSegDist(const stfPoint_t& p1,const stfPoint_t& p2,const stfPoint_t& p);
	double GetPointDistance(const stfPoint_t& p1, const stfPoint_t& p2);
        bool  doubleEqual(double a, double b);
protected:
	unsigned int fPointID;
	unsigned int fEdgeID;

};

#endif CTRAFICGRAPH_H_
