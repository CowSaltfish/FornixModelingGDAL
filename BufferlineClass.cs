using System;

namespace FornixModelingGDAL
{
    /// <summary>
    /// 计算水系的缓冲区边界
    /// </summary>
    class BufferlineClass
    {
        FeatureDispose fd = new FeatureDispose();
        InputData id = new InputData();
        DataTrans dt = new DataTrans();
        /// <summary>
        /// 执行线缓冲区分析
        /// </summary>
        /// <param name="waterlineFile">输入水系shp文件</param>
        /// <param name="outfolder">选择输出文件位置</param>
        /// <param name="radius">缓冲区半径</param>
        public void BufferlineAnalyse(string waterlineFile, string outfolder, double radius)
        {
            id.DeleteFiles(outfolder);
            IFeatureClass featureClass = id.InputShp(waterlineFile);//读取水系数据
            for (int i = 0; i < featureClass.FeatureCount(null); i++)//遍历水系要素
            {
                IPolyline waterline = dt.FeaturetoPolyline(featureClass, i);//获取一条水系要素,将水系矢量数据转换为图形数据
                IPointCollection pointCollection = waterline as IPointCollection;//将水系矢量数据转换为点集图形数据  
                GetBufferEdgeCoords(outfolder, pointCollection, radius, i);
            }
        }
        /// <summary>
        /// 根据给定的一系列有顺序的坐标，逆时针生成缓冲区的边界坐标。
        /// </summary>
        /// <param name="strPolyLineCoords">一系列有顺序的坐标</param>
        /// <param name="radius">缓冲区半径</param>
        /// <returns>缓冲区的边界坐标</returns>
        public void GetBufferEdgeCoords(VertexCollection polyline, double radius, int i)
        {
            if (polyline.Count > 0)
            {
                //分别生成内环和外环的缓冲区边界点坐标串
                VertexCollection lpc = GetLeftBufferEdgeCoords(polyline, radius);
                VertexCollection polyliner = Reverse(polyline);
                VertexCollection rpc = GetLeftBufferEdgeCoords(polyliner, radius);
                VertexCollection rpcr = Reverse(rpc);
            }
        }
        /// <summary>
        /// 执行线要素的反转
        /// </summary>
        /// <param name="ipc"></param>
        /// <returns></returns>
        public VertexCollection Reverse(VertexCollection ipc)
        {
            VertexCollection ipc0 = new VertexCollection();
            for (int i = ipc.Count - 1; i >= 0; i--)
            {
                ipc0.addVer(ipc.getVer(i));
            }
            return ipc0;
        }
        /// <summary>
        /// 计算向量方位角
        /// </summary>
        /// <param name="preCoord">向量起点</param>
        /// <param name="nextCoord">向量终点</param>
        /// <returns>返回弧度角</returns>
        public double GetQuadrantAngle(Vertex preCoord, Vertex nextCoord)
        {
            return GetQuadrantAngle(nextCoord.X() - preCoord.X(), nextCoord.Y() - preCoord.Y());
        }
        /// <summary>
        /// 由增量X和增量Y所形成的向量的象限角度
        /// </summary>
        /// <param name="x">增量X</param>
        /// <param name="y">增量Y</param>
        /// <returns>象限角</returns>
        public double GetQuadrantAngle(double x, double y)
        {
            double theta = Math.Atan(y / x);
            if (x > 0 && y < 0) return Math.PI * 2 + theta;
            if (x < 0) return theta + Math.PI;
            return theta;
        }
        /// <summary>
        /// 获取由相邻的三个点所形成的两个向量之间的夹角
        /// </summary>
        /// <param name="preCoord"></param>
        /// <param name="midCoord"></param>
        /// <param name="nextCoord"></param>
        /// <returns></returns>
        public double GetIncludedAngle(Vertex preCoord, Vertex midCoord, Vertex nextCoord)
        {
            double innerProduct = (midCoord.X() - preCoord.X()) * (nextCoord.X() - midCoord.X()) + (midCoord.Y() - preCoord.Y()) * (nextCoord.Y() - midCoord.Y());
            double mode1 = Math.Sqrt(Math.Pow((midCoord.X() - preCoord.X()), 2.0) + Math.Pow((midCoord.Y() - preCoord.Y()), 2.0));
            double mode2 = Math.Sqrt(Math.Pow((nextCoord.X() - midCoord.X()), 2.0) + Math.Pow((nextCoord.Y() - midCoord.Y()), 2.0));
            double acos = innerProduct / (mode1 * mode2);
            if (acos > 1) acos = 1;
            if (acos < -1) acos = -1;
            return Math.Acos(acos);
        }

        /// <summary>
        /// 根据给定的一系列有顺序的坐标，逆时针生成轴线左侧的缓冲区边界点
        /// </summary>
        /// <param name="coords">一系列有顺序的坐标</param>
        /// <param name="radius">缓冲区半径</param>
        /// <returns>缓冲区的边界坐标</returns>
        public VertexCollection GetLeftBufferEdgeCoords(VertexCollection coords, double radius)
        {
            VertexCollection polyline = new VertexCollection();
            Vertex point = new Vertex();
            //计算时所需变量
            double alpha = 0.0;//向量绕起始点沿顺时针方向旋转到X轴正半轴所扫过的角度
            double delta = 0.0;//前后线段所形成的向量之间的夹角
            double l = 0.0;//前后线段所形成的向量的叉积
            //辅助变量
            double startRadian = 0.0;
            double endRadian = 0.0;
            double beta = 0.0;
            double x = 0.0, y = 0.0;
            //中间节点
            for (int i = 1; i < coords.Count - 1; i++)
            {
                alpha = GetQuadrantAngle(coords.getVer(i), coords.getVer(i + 1));
                delta = GetIncludedAngle(coords.getVer(i - 1), coords.getVer(i), coords.getVer(i + 1));
                l = GetVectorProduct(coords.getVer(i - 1), coords.getVer(i), coords.getVer(i + 1));
                if (l > 0)//凸
                {
                    startRadian = alpha + (3 * Math.PI) / 2 - delta;
                    endRadian = alpha + (3 * Math.PI) / 2;
                    VertexCollection ipc1 = GetBufferCoordsByRadian(coords.getVer(i), startRadian, endRadian, radius);
                    for (int j = 0; j < ipc1.Count; j++)
                    {
                        if (ipc1.getVer(j).X() > 0)
                        { }
                        Alter(ref polyline, ipc1.getVer(j));
                        polyline.addVer(ipc1.getVer(j));
                    }
                }
                else if (l < 0)
                {
                    beta = alpha - (Math.PI - delta) / 2;
                    x = Math.Round(coords.getVer(i).X() + radius * Math.Cos(beta), 2);
                    y = Math.Round(coords.getVer(i).Y() + radius * Math.Sin(beta), 2);
                    Vertex ipoint = new Vertex();
                    ipoint.X(x);
                    ipoint.Y(y);
                    if (ipoint.X() > 0) { }
                    Alter(ref polyline, ipoint);
                    polyline.addVer(ipoint);
                }
            }
            return polyline;
        }
        /// <summary>
        /// 判断新加入的点会不会是边界产生自相交，是则进行修改
        /// </summary>
        /// <param name="ipc"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        public bool Alter(ref VertexCollection ipc, Vertex point)
        {
            int count = ipc.Count;
            if (count >= 3)
            {
                bool flag = false;
                int index = count - 1;
                for (int i = index - 1; i >= 1; i--)
                {

                    if (InsectionJudge(ipc.getVer(index), point, ipc.getVer(i), ipc.getVer(i - 1)))
                    {
                        flag = true;
                        index = i;
                        break;
                    }
                }
                if (flag)
                {
                    Vertex insectp = Inter(ipc.getVer(count - 1), point, ipc.getVer(index), ipc.getVer(index - 1));
                    ipc.removeVers(index, count - index);
                    ipc.addVer(insectp);
                }
            }
            return true;
        }

        /// <summary>
        /// 判断俩直线是否相交
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        public bool InsectionJudge(Vertex a, Vertex b, Vertex c, Vertex d)
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
        /// 计算交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        /// <returns></returns>
        public Vertex Inter(Vertex p1, Vertex p2, Vertex p3, Vertex p4)
        {
            Vertex point = new Vertex();
            double s1 = fArea(p1, p2, p3), s2 = fArea(p1, p2, p4);
            point.X(Math.Round((p4.X() * s1 + p3.X() * s2) / (s1 + s2), 2));
            point.Y(Math.Round((p4.Y() * s1 + p3.Y() * s2) / (s1 + s2), 2));
            return point;
        }

        public double Cross(Vertex p1, Vertex p2, Vertex p3, Vertex p4)
        {
            return (p2.X() - p1.X()) * (p4.Y() - p3.Y()) - (p2.Y() - p1.Y()) * (p4.X() - p3.X());
        }

        public double Area(Vertex p1, Vertex p2, Vertex p3)
        {
            return Cross(p1, p2, p1, p3);
        }

        public double fArea(Vertex p1, Vertex p2, Vertex p3)
        {
            return Math.Abs(Area(p1, p2, p3));
        }
        /// <summary>
        /// 获取指定弧度范围之间的缓冲区圆弧拟合边界点
        /// </summary>
        /// <param name="center">指定拟合圆弧的原点</param>
        /// <param name="startRadian">开始弧度</param>
        /// <param name="endRadian">结束弧度</param>
        /// <param name="radius">缓冲区半径</param>
        /// <returns>缓冲区的边界坐标</returns>
        private VertexCollection GetBufferCoordsByRadian(Vertex center, double startRadian, double endRadian, double radius)
        {
            VertexCollection points = new VertexCollection();
            double gamma = Math.PI / 100;
            double x = 0.0, y = 0.0;
            for (double phi = startRadian; phi <= endRadian + 0.000000000000001; phi += gamma)
            {
                Vertex point = new Vertex();
                x = Math.Round(center.X() + radius * Math.Cos(phi), 2);
                y = Math.Round(center.Y() + radius * Math.Sin(phi), 2);
                point.X(x);
                point.Y(y);
                points.addVer(point);
            }
            return points;
        }
        /// <summary>
        /// 获取相邻三个点所形成的两个向量的交叉乘积
        /// </summary>
        /// <param name="preCoord">第一个节点坐标</param>
        /// <param name="midCoord">第二个节点坐标</param>
        /// <param name="nextCoord">第三个节点坐标</param>
        /// <returns>相邻三个点所形成的两个向量的交叉乘积</returns>
        private double GetVectorProduct(Vertex preCoord, Vertex midCoord, Vertex nextCoord)
        {
            return (midCoord.X() - preCoord.X()) * (nextCoord.Y() - midCoord.Y()) - (nextCoord.X() - midCoord.X()) * (midCoord.Y() - preCoord.Y());
        }
    }
}
