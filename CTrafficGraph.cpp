#include "CTrafficGraph.h"
#include <cmath>
bool doubleEqual(double a, double b)
{
	if (fabs(a - b) < 0.0001) {
		return true;
	}

	return false;
}

double GetPointDistance(const stfPoint_t& p1, const stfPoint_t& p2)
{
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

//计算点到线段距离
double PointToSegDist(const stfPoint_t& p1, const stfPoint_t& p2, const stfPoint_t& p)
{
	float a, b, c;
	a = GetPointDistance(p2, p);
	if (a <= 0.00001)
		return 0.0f;
	b = GetPointDistance(p1, p);
	if (b <= 0.00001)
		return 0.0f;
	c = GetPointDistance(p1, p2);
	if (c <= 0.00001)
		return a;//如果p1和p2坐标相同，则退出函数，并返回距离 

		/*钝角 Obtuse angle
			
		  p2*
       \ *
				\  * 
				 \____*p
				 p1
				 */
	float aa = a * a;
	float bb = b * b;
	float cc = c * c;
	if (a*a > b * b + c * c + 0.00001)
		return 0.0f;
		/*钝角 Obtuse angle
	   p1*
       \ *
				\  * 
				 \____*p
				 p2
	  */
	if (b*b  > a * a + c * c + 0.00001)
		return 0.0f;

	float l = (a + b + c) / 2;//周长的一半
	float s = sqrt(l*(l - a)*(l - b)*(l - c));//海伦公式求面积
	return 2 * s / c;
}

CTrafficGraph::CTrafficGraph(vector<stfEdge>& vecEdges)
{
	if (vecEdges.empty())
	{
		return;
	}
	fPointID = 0;
	fEdgeID = 0;
	
	for (int i = 0; i < vecEdges.size(); i++)
	{
		InputEdge(vecEdges[i].StartPoint, vecEdges[i].EndPoint);
	}
}

CTrafficGraph::~CTrafficGraph()
{
}

bool CTrafficGraph::hasPoint(const stfPoint_t& Point, unsigned int& ID)
{
	map < stfPoint_t, unsigned int >::iterator it = fmapPoints.begin();
	for (; it != fmapPoints.end(); it++) {
		if (it->first == Point) {
			ID = it->second;
			return true;
		}
	}
	return false;
}

bool CTrafficGraph::hasEdge(const stfPoint_t& startPoint,const stfPoint_t& endPoint, unsigned int& ID)
{
	stfEdge_t tmpEdge;
	tmpEdge.StartPoint = startPoint;
	tmpEdge.EndPoint = endPoint;
	map < stfEdge_t, unsigned int >::iterator it = fmapEdges.begin();
	for (; it != fmapEdges.end(); it++) {
		if (it->first == tmpEdge) {
			ID = it->second;
			return true;
		}
	}
	return false;
}

bool CTrafficGraph::PointEqual(const stfPoint_t& Point1,const stfPoint_t& Point2)
{
	if (Point1 == Point2) {
		return true;
	}

	return false;
}

bool IsInLineSegment(stfPoint_t& crossPoint, stfPoint_t& startPoint, stfPoint_t& endPoint)
{
	double x1, y1, x2, y2;
	if (startPoint.x > endPoint.x) {
		x1 = endPoint.x;
		x2 = startPoint.x;
	}
	else {
		x1 = startPoint.x;
		x2 = endPoint.x;
	}

	if (startPoint.y > endPoint.y) {
		y1 = endPoint.y;
		y2 = startPoint.y;
	}
	else {
		y1 = startPoint.y;
		y2 = endPoint.y;
	}
	bool flag1 = false;
	bool flag2 = false;
	//判断是否是垂直线
	if (doubleEqual(x1,x2))
	{
		if (doubleEqual(x1,crossPoint.x))
		{
			flag1 = true;
		}
		else
			return false;
	}
	else
	{
		if ((crossPoint.x >= x1) && (crossPoint.x <= x2))
		{
			flag1 = true;
		}
		else
			return false;

	}

	if (doubleEqual(y1,y2))
	{
		if (doubleEqual(y1,crossPoint.y))
		{
			flag2 = true;
		}
		else 
			return false;
	}
	else
	{
		if ( (crossPoint.y >= y1) && (crossPoint.y <= y2))
		{
			flag2  = true;
		}
		else
			return false;
	}
	if (flag1 && flag2)
		return true;

	return false;
}

// 判断两条线段是否相交，如果相交，返回交点
bool CTrafficGraph::lineIntersection(stfPoint_t& startPoint, stfPoint_t& endPoint, stfEdge_t& edge, stfPoint_t& crossPoint)
{
	double x1 = startPoint.x;
	double y1 = startPoint.y;
	double x2 = endPoint.x;
	double y2 = endPoint.y;

	double x3 = edge.StartPoint.x;
	double y3 = edge.StartPoint.y;
	double x4 = edge.EndPoint.x;
	double y4 = edge.EndPoint.y;

	if (PointEqual(startPoint, edge.StartPoint))
		return false;
	if (PointEqual(startPoint, edge.EndPoint))
		return false;
	if (PointEqual(endPoint, edge.StartPoint))
		return false;
	if (PointEqual(endPoint, edge.EndPoint))
		return false;

	// 判断是否是竖直线（特殊斜率情况）
	bool gt1 = doubleEqual(x2 - x1, 0);
	bool gt2 = doubleEqual(x4 - x3, 0);
	if (gt1 && gt2)
		return false;

	double gradient1 = 0;
	double gradient2 = 0;

	if (!gt1)
		gradient1 = (y2 - y1) / (x2 - x1);
	else {
		gradient2 = (y4 - y3) / (x4 - x3);
		double a2 = gradient2;
		double b2 = y3 - a2 * x3;

		double xc = x1;
		double yc = a2 * xc + b2;

		crossPoint.x = xc;
		crossPoint.y = yc;

		bool bin1 = IsInLineSegment(crossPoint, startPoint, endPoint);
		bool bin2 = IsInLineSegment(crossPoint, edge.StartPoint, edge.EndPoint);
		if (bin1 && bin2)
			return true;
		else
			return false;
	}

	if (!gt2)
		gradient2 = (y4 - y3) / (x4 - x3);
	else {
		gradient1 = (y2 - y1) / (x2 - x1);
		double a1 = gradient1;
		double b1 = y1 - a1 * x1;
		double xc = x3;
		double yc = a1 * xc + b1;

		crossPoint.x = xc;
		crossPoint.y = yc;

		bool bin1 = IsInLineSegment(crossPoint, startPoint, endPoint);
		bool bin2 = IsInLineSegment(crossPoint, edge.StartPoint, edge.EndPoint);
		if (bin1 && bin2)
			return true;
		else
			return false;
	}

	// 如果两根线都不是竖直线，判断是否平行
	if (doubleEqual(gradient1, gradient2))
		return false;

	double a1 = gradient1;
	double b1 = y1 - a1 * x1;

	double a2 = gradient2;
	double b2 = y3 - a2 * x3;

	double xc = (b2 - b1) / (a1 - a2);
	double yc = a1 * xc + b1;

	crossPoint.x = xc;
	crossPoint.y = yc;

	bool bin1 = IsInLineSegment(crossPoint, startPoint, endPoint);
	bool bin2 = IsInLineSegment(crossPoint, edge.StartPoint, edge.EndPoint);
	if (bin1 && bin2)
		return true;
	else
		return false;
}

unsigned int CTrafficGraph::AddPoint(stfPoint_t& Point)
{
	unsigned int ID = 0;
	if (!hasPoint(Point, ID)) {
		ID = ++fPointID;
		Point.pid = ID;
		fmapPoints[Point] = ID;
	}
	else
	{
		Point.pid = ID;
	}
	return ID;
}

unsigned int  CTrafficGraph::AddEdge( stfPoint_t& startPoint, stfPoint_t& endPoint)
{
	unsigned int startPID = 0;
	unsigned int endPID = 0;
	startPID = AddPoint(startPoint);
	endPID = AddPoint(endPoint);

	stfEdge_t tmpEdge;
	tmpEdge.StartPoint = startPoint;
	tmpEdge.EndPoint = endPoint;
	
	unsigned int EdgeID = 0;
	if (hasEdge(startPoint, endPoint, EdgeID))
		return EdgeID;

	tmpEdge.eid = ++fEdgeID;
	tmpEdge.Weight = 0;
	fmapEdges[tmpEdge] = tmpEdge.eid;
	return tmpEdge.eid;
	
}

void CTrafficGraph::InputEdge(stfPoint_t& startPoint, stfPoint_t& endPoint)
{
	unsigned int startPID = 0;
	unsigned int EndPID = 0;
	bool retStart = hasPoint(startPoint, startPID);
	bool retEnd = hasPoint(endPoint, EndPID);

	//若两点为同一个点，直接退出
	if (startPoint == endPoint)
		return;

	//若两点已经存在，直接退出
	/*if (retStart && retEnd)
		return;*/

	unsigned int EdgeID = 0;
	if (hasEdge(startPoint, endPoint, EdgeID))
		return;

	// 边信息
	bool hasCross = false;
	stfPoint_t crossPoint;

	map < stfEdge_t, unsigned int >::iterator it = fmapEdges.begin();
	for (; it != fmapEdges.end(); it++) {
		stfEdge_t tet = it->first;
		if (lineIntersection(startPoint, endPoint, tet, crossPoint)) {
			unsigned int pID = AddPoint(crossPoint);
			hasCross = true;

			// 如果交点为已有线段的端点，则已有线段无需打断
			if ((crossPoint == tet.StartPoint) || (crossPoint == tet.EndPoint)) {
				break;
			}

			// 添加打断的老线段后生成的两个线段
			AddEdge(tet.StartPoint, crossPoint);
			AddEdge(crossPoint, tet.EndPoint);
			// 删除原来的老线段
			fmapEdges.erase(it);
			break;
		}
	}

	if (hasCross) {
		// 新加入的边被crosspoint打断生成的子线段继续InputEdge
		InputEdge(startPoint, crossPoint);
		InputEdge(crossPoint, endPoint);
	}
	else {
		AddEdge(startPoint, endPoint);
	}
	//DumpEdge();
}

void CTrafficGraph::DumpEdge()
{
	map < stfEdge_t, unsigned int >::iterator it = fmapEdges.begin();
	for (; it != fmapEdges.end(); it++) {
		stfPoint_t p1 = it->first.StartPoint;
		stfPoint_t p2 = it->first.EndPoint;
		double dis = sqrt((p2.y - p1.y)*(p2.y - p1.y) + (p2.x - p1.x)*(p2.x - p1.x));
		//printf("Edge:  point1.x:%.6f point1.y:%.6f  point2.x:%.6f point2.y:%.6f---PID1:%d,PID2:%d-- EID:%d--- distance:%.6f\n", p1.x, p1.y, p2.x, p2.y,p1.pid,p2.pid,it->second,dis);
		printf("Edge:  point1:(%.3f, %.3f)  point2:(%.3f, :%.3f)\n", p1.x, p1.y, p2.x, p2.y);
	}
}
