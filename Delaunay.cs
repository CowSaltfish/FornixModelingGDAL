using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace FornixModelingGDAL
{
    /// <summary>
    /// 专门用于生成2D穹窿状Delaunay三角网
    /// 先插入外环顶点生成初始三角网
    /// 再插入内环顶点
    /// 对三角网进行边界恢复与调整（法向调整、删除多余的三角面），生成内外环控制的三角网
    /// 根据包含关系，插入内部点
    /// 再一次调整法向
    /// </summary>
    class Delaunay
    {
        /// <summary>
        /// 内部点序号
        /// </summary>
        private static int _idxInnerPoint = 0;
        /// <summary>
        /// 顶面或底面
        /// </summary>
        private static bool _UpOrDown = false;
        /// <summary>
        /// 三角网
        /// </summary>
        private static SortedList<string, Triangle> _tris;
        /// <summary>
        /// 空腔
        /// </summary>
        private static SortedList<string, Triangle> _cavity;
        /// <summary>
        /// 空腔边界
        /// </summary>
        private static List<Edge> _edges;
        /// <summary>
        /// 超级矩形：左上，右上，左下，右下
        /// </summary>
        private static Vertex[] _rectangle = new Vertex[4];

        /// <summary>
        /// 根据多边形内外环、内部点生成Delauney三角网
        /// </summary>
        /// <param name="tris"></param>
        /// <param name="outRing"></param>
        /// <param name="inRing"></param>
        /// <param name="inner"></param>
        /// <returns>多边形内部点点集</returns>
        public static VertexCollection GenTriMesh(SortedList<string, Triangle> tris, VertexCollection outRing, VertexCollection inRing, VertexCollection inner)
        {
            _tris = tris;
            VertexCollection vcInner;
            if (null != inner) _UpOrDown = true;
            else _UpOrDown = false;

            //建立矩形包围点集
            CreateSuperRectangle(outRing);
            //将多边形外环顶点插入三角网中
            InsertVC2Tris(outRing);
            //将多边形内环顶点插入三角网中
            InsertVC2Tris(inRing);
            //删去多余三角形，调整法向
            PostTreatment(outRing, inRing);
            //边界恢复
            //EdgesRecover(outRing, inRing);
            //将多边形内部顶点插入三角网中
            vcInner = InsertVC2Tris(inner);
            return vcInner;
        }

        /// <summary>
        /// 建立矩形包围点集
        /// </summary>
        /// <param name="outRing"></param>
        private static void CreateSuperRectangle(VertexCollection outRing)
        {
            _rectangle[0] = new Vertex();
            _rectangle[1] = new Vertex();
            _rectangle[2] = new Vertex();
            _rectangle[3] = new Vertex();
            double minX = 100000000, maxX = 0.0, minY = 100000000, maxY = 0.0;
            for (int i = 0; i < outRing.Count; ++i)
            {
                Vertex ver = outRing.getVer(i);
                if (ver.X() < minX) minX = ver.X();
                if (ver.X() > maxX) maxX = ver.X();
                if (ver.Y() < minY) minY = ver.Y();
                if (ver.Y() > maxY) maxY = ver.Y();
            }
            _rectangle[0].X(minX - (maxX - minX) / 10); _rectangle[0].Y(maxY + (maxY - minY) / 10); _rectangle[0].ID = -1;
            _rectangle[1].X(maxX + (maxX - minX) / 10); _rectangle[1].Y(maxY + (maxY - minY) / 10); _rectangle[1].ID = -2;
            _rectangle[2].X(minX - (maxX - minX) / 10); _rectangle[2].Y(minY - (maxY - minY) / 10); _rectangle[2].ID = -3;
            _rectangle[3].X(maxX + (maxX - minX) / 10); _rectangle[3].Y(minY - (maxY - minY) / 10); _rectangle[3].ID = -4;
        }

        /// <summary>
        /// 由多边形顶点生成三角网
        /// </summary>
        /// <param name="vc"></param>
        private static VertexCollection InsertVC2Tris(VertexCollection vc)
        {
            if (null == vc)
                return null;
            Triangle tri;
            _cavity = new SortedList<string, Triangle>();
            _edges = new List<Edge>();
            int i = 0, j;
            VertexCollection vcInner = new VertexCollection();
            if (0 == _tris.Count)
            {
                //第一个点连接矩形顶点
                Vertex ver0 = vc.getVer(0);
                Triangle[] triOfRectangle = new Triangle[4];
                triOfRectangle[0] = new Triangle();
                triOfRectangle[1] = new Triangle();
                triOfRectangle[2] = new Triangle();
                triOfRectangle[3] = new Triangle();
                triOfRectangle[0].addPoint(ver0);
                triOfRectangle[0].addPoint(_rectangle[1]);
                triOfRectangle[0].addPoint(_rectangle[0]);
                _tris.Add(triOfRectangle[0].name, triOfRectangle[0]);
                triOfRectangle[1].addPoint(ver0);
                triOfRectangle[1].addPoint(_rectangle[3]);
                triOfRectangle[1].addPoint(_rectangle[1]);
                _tris.Add(triOfRectangle[1].name, triOfRectangle[1]);
                triOfRectangle[2].addPoint(ver0);
                triOfRectangle[2].addPoint(_rectangle[2]);
                triOfRectangle[2].addPoint(_rectangle[3]);
                _tris.Add(triOfRectangle[2].name, triOfRectangle[2]);
                triOfRectangle[3].addPoint(ver0);
                triOfRectangle[3].addPoint(_rectangle[0]);
                triOfRectangle[3].addPoint(_rectangle[2]);
                _tris.Add(triOfRectangle[3].name, triOfRectangle[3]);
                ++i;
            }
            //将地层顶点插入三角网中
            for (; i < vc.Count; ++i)
            {
                Vertex verInserted = vc.getVer(i);
                if (verInserted.toRemove)
                    continue;
                if (!verInserted.innerPoint)
                {
                    for (j = 0; j < _tris.Count; ++j)
                    {
                        tri = _tris.Values[j];
                        Vertex circumcenter = tri.calCircumCenter();
                        //点在三角形外接圆中
                        if (verInserted.calDistance(circumcenter) <= tri.points.getVer(0).calDistance(circumcenter))
                        {
                            //广度遍历得到空腔
                            FindCavity(verInserted, tri);
                            //由空腔生成新三角形
                            CreateTriByCavity(verInserted);
                            break;
                        }
                    }
                }
                else
                {
                    for (j = 0; j < _tris.Count; ++j)
                    {
                        tri = _tris.Values[j];
                        //判断点在三角形中
                        if (verInserted.inside(tri))
                        {
                            verInserted.ID = ++_idxInnerPoint;
                            vcInner.addVer(verInserted);
                            //广度遍历得到空腔
                            FindCavity(verInserted, tri);
                            //由空腔生成新三角形
                            CreateTriByCavity(verInserted);
                            break;
                        }
                    }
                }

            }
            //插入内部点后调整三角面法向
            if (0 != vc.Count && vc.getVer(0).innerPoint)
            {
                for (i = 0; i < _tris.Count; ++i)
                {
                    tri = _tris.Values[i];
                    if (tri.calVector() < 0.0)
                        tri.reverse();
                }
            }
            return vcInner;
        }

        /// <summary>
        /// 后处理，删去多余三角形，调整法向
        /// </summary>
        /// <param name="outRing"></param>
        /// <param name="inRing"></param>
        private static void PostTreatment(VertexCollection outRing, VertexCollection inRing)
        {
            int minOutRingID = outRing.getVer(0).ID;
            int maxOutRingID = outRing.getVer(outRing.Count - 1).ID;
            int minInRingID = 0 != inRing.Count ? inRing.getVer(0).ID : 0;
            int maxInRingID = 0 != inRing.Count ? inRing.getVer(inRing.Count - 1).ID : 0;
            Triangle tri;
            int i;
            List<String> nameOFTriRemoved = new List<string>();

            for (i = 0; i < _tris.Count; ++i)
            {
                tri = _tris.Values[i];
                //删去与超级矩形顶点有关的三角形、内环内的三角形、外环外的三角形
                if (tri.hasVer(-1) || tri.hasVer(-2) || tri.hasVer(-3) || tri.hasVer(-4)
                    || (0 != inRing.Count && tri.pointsOnRing(minInRingID, maxInRingID) && tri.calVector() < 0.0)
                    || (tri.pointsOnRing(minOutRingID, maxOutRingID) && tri.calVector() > 0.0))
                    nameOFTriRemoved.Add(tri.name);
                //调整三角形法向
                else if (!_UpOrDown && tri.calVector() > 0.0)
                    tri.reverse();
                else if (_UpOrDown && tri.calVector() < 0.0)
                    tri.reverse();
            }
            //执行删除
            for (i = 0; i < nameOFTriRemoved.Count; ++i)
                _tris.Remove(nameOFTriRemoved[i]);
        }

        /// <summary>
        /// 边界恢复
        /// </summary>
        /// <param name="outRing"></param>
        /// <param name="inRing"></param>
        private static void EdgesRecover(VertexCollection outRing, VertexCollection inRing)
        {
            Vertex ver0, ver1;
            Triangle tri;
            List<Triangle> trisWithVer0 = new List<Triangle>();
            int i, j;
            bool flag = false;
            for (i = 0; i < outRing.Count - 1; ++i)
            {
                ver0 = outRing.getVer(i);
                ver1 = outRing.getVer(i + 1);
                trisWithVer0.Clear();
                for (j = 0; j < _tris.Count; ++j)
                {
                    flag = false;
                    tri = _tris.Values[j];
                    if (tri.hasVer(ver0) && tri.hasVer(ver1))
                    {
                        flag = true;
                        break;
                    }
                    else if (tri.hasVer(ver0))
                    {
                        trisWithVer0.Add(tri);
                    }
                }
                if (!flag)//在三角网中没有找到这个边
                {
                    for (j = 0; j < trisWithVer0.Count; ++j)
                    {
                        tri = trisWithVer0[j];
                        //目标线段与该三角形第三条边相交
                        if (InsectionJudge(ver0, ver1, tri.points.getVer(1), tri.points.getVer(1)))
                        {
                            //广度遍历该三角形，找到所有与目标线段相交的三角形，形成一个空腔
                            FindERCavity(tri, ver0, ver1);
                            //普里姆算法，先连接目标线段，再依次将离目标线段（多边形）最近的顶点加入三角网中，直至填满空腔
                            RecoverTris();
                        }
                    }
                }
            }
        }


        /// <summary>
        /// 寻找空腔
        /// </summary>
        /// <param name="verInserted"></param>
        /// <param name="tri"></param>
        private static void FindCavity(Vertex verInserted, Triangle tri)
        {
            Queue<Triangle> triQueue = new Queue<Triangle>();
            _cavity.Clear();
            _edges.Clear();
            triQueue.Clear();
            triQueue.Enqueue(tri);

            //广度遍历生成空腔
            while (0 != triQueue.Count)
            {
                Triangle curTri = triQueue.Dequeue();
                if (!curTri.isCavity)
                    AddTri2Cavity(curTri);
                foreach (Triangle adjtri in curTri.getAdjTri(_tris))
                {
                    if (adjtri.isCavity)
                        continue;
                    Vertex circumcenter = adjtri.calCircumCenter();
                    if (verInserted.calDistance(circumcenter) <= adjtri.points.getVer(0).calDistance(circumcenter))
                        triQueue.Enqueue(adjtri);
                }
            }
        }

        /// <summary>
        /// 将三角形加入空腔
        /// </summary>
        /// <param name="curTri"></param>
        private static void AddTri2Cavity(Triangle curTri)
        {
            Edge edge, newEdge0, newEdge1, newEdge2;
            newEdge0 = new Edge(curTri.points.getVer(0), curTri.points.getVer(1));
            newEdge1 = new Edge(curTri.points.getVer(1), curTri.points.getVer(2));
            newEdge2 = new Edge(curTri.points.getVer(0), curTri.points.getVer(2));
            List<Edge> edges1 = new List<Edge>();
            bool flag;

            for (int i = 0; i < _edges.Count; ++i)
            {
                edge = _edges[i];
                flag = false;
                if (edge.Equals(curTri.points.getVer(0), curTri.points.getVer(1)))
                {
                    flag = true;
                    newEdge0 = null;
                }
                else if (edge.Equals(curTri.points.getVer(1), curTri.points.getVer(2)))
                {
                    flag = true;
                    newEdge1 = null;
                }
                else if (edge.Equals(curTri.points.getVer(2), curTri.points.getVer(0)))
                {
                    flag = true;
                    newEdge2 = null;
                }
                if (!flag)
                    edges1.Add(edge);
            }
            if (null != newEdge0)
                edges1.Add(newEdge0);
            if (null != newEdge1)
                edges1.Add(newEdge1);
            if (null != newEdge2)
                edges1.Add(newEdge2);
            _cavity.Add(curTri.name, curTri);
            curTri.isCavity = true;
            _edges.Clear();
            _edges.AddRange(edges1);
        }

        /// <summary>
        /// 在空腔中新建三角面
        /// </summary>
        /// <param name="verInserted"></param>
        private static void CreateTriByCavity(Vertex verInserted)
        {
            Triangle newTri;
            int i;
            string nameR;
            //删去空腔中三角形
            for (i = 0; i < _cavity.Count; ++i)
            {
                if (_tris.ContainsKey(_cavity.Keys[i]))
                    _tris.Remove(_cavity.Keys[i]);
                else
                {
                    nameR = _cavity.Values[i].reverseName();
                    _tris.Remove(nameR);
                }
            }
            //添加新三角形
            for (i = 0; i < _edges.Count; ++i)
            {
                newTri = new Triangle();
                newTri.addPoint(_edges[i].v0);
                newTri.addPoint(_edges[i].v1);
                newTri.addPoint(verInserted);
                _tris.Add(newTri.name, newTri);
            }
        }

        /// <summary>
        /// 边界恢复中寻找空腔
        /// </summary>
        /// <param name="tri"></param>
        /// <param name="ver0"></param>
        /// <param name="ver1"></param>
        private static void FindERCavity(Triangle tri, Vertex ver0, Vertex ver1)
        {
            Queue<Triangle> triQueue = new Queue<Triangle>();
            _cavity.Clear();
            _edges.Clear();
            triQueue.Clear();
            triQueue.Enqueue(tri);

            //广度遍历生成空腔
            while (0 != triQueue.Count)
            {
                Triangle curTri = triQueue.Dequeue();
                if (!curTri.isCavity)
                    AddTri2Cavity(curTri);
                foreach (Triangle adjtri in curTri.getAdjTri(_tris))
                {
                    if (adjtri.isCavity)
                        continue;
                    //三角形包含ver1
                    else if(true)
                        ;
                    //三角形被目标线段穿过
                    else if (true)
                        triQueue.Enqueue(adjtri);
                }
            }
        }


        private static void RecoverTris()
        {
            
        }

        /// <summary>
        /// 判断俩直线是否相交
        /// </summary>
        /// <param name="a">第一条直线第一个顶点</param>
        /// <param name="b">第一条直线第二个顶点</param>
        /// <param name="c">第二条直线第一个顶点</param>
        /// <param name="d">第二条直线第二个顶点</param>
        /// <returns></returns>
        private static bool InsectionJudge(Vertex a, Vertex b, Vertex c, Vertex d)
        {
            
            /*
            快速排斥：
            两个线段为对角线组成的矩形，如果这两个矩形没有重叠的部分，那么两条线段是不可能出现重叠的
            */
            if (!(Math.Min(a.X(), b.X()) <= Math.Max(c.X(), d.X()) &&
                Math.Min(c.Y(), d.Y()) <= Math.Max(a.Y(), b.Y()) &&
                Math.Min(c.X(), d.X()) <= Math.Max(a.X(), b.X()) &&
                Math.Min(a.Y(), b.Y()) <= Math.Max(c.Y(), d.Y())))//这一步是判定两矩形是否相交
            {
                return false;
            }
            /*
            跨立实验：
            如果两条线段相交，那么必须跨立，就是以一条线段为标准，另一条线段的两端点一定在这条线段的两段
            也就是说a b两点在线段cd的两端，c d两点在线段ab的两端
            */
            double u, v, w, z;//分别记录两个向量
            u = (c.X() - a.X()) * (b.Y() - a.Y()) - (b.X() - a.X()) * (c.Y() - a.Y());
            v = (d.X() - a.X()) * (b.Y() - a.Y()) - (b.X() - a.X()) * (d.Y() - a.Y());
            w = (a.X() - c.X()) * (d.Y() - c.Y()) - (d.X() - c.X()) * (a.Y() - c.Y());
            z = (b.X() - c.X()) * (d.Y() - c.Y()) - (d.X() - c.X()) * (b.Y() - c.Y());
            return (u * v <= 0.00000001 && w * z <= 0.00000001);
        }

        /// <summary>
        /// 求线段与当前线段的交点
        /// </summary>
        /// <param name="edge">目标线段</param>
        /// <returns>交点</returns>
        private static Vertex GetCrossPoint(Vertex a, Vertex b, Vertex c, Vertex d)
        {
            Vertex crossPoint = new Vertex();
            double A = (b.X() - a.X()) / (b.Y() - a.Y());
            double B = (d.X() - c.X()) / (d.Y() - c.Y());
            double y = (c.X() - a.X() + A * a.Y() - B * c.Y()) / (A - B);
            double x = A * y - A * a.Y() + a.X();
            crossPoint.X(x);
            crossPoint.Y(y);
            crossPoint.Z((a.Z() + b.Z() + c.Z() + d.Z()) / 4);
            crossPoint.ID = b.ID;
            return crossPoint;
        }
    }
}
