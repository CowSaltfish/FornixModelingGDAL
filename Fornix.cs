using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using OSGeo.GDAL;

namespace FornixModelingGDAL
{
    class Fornix
    {
        public string name { get; set; }
        public UpFace upFace { get; private set; }
        public DownFace downFace { get; private set; }
        public OutSide outSide { get; private set; }
        public InSide inSide { get; private set; }
        public double dip { get; set; }
        public double[] range { get; private set; }
        public int Bessel_num_collected;//贝塞尔采集点数

        public Fornix()
        {
            this.upFace = new UpFace();
            this.downFace = new DownFace();
            this.outSide = new OutSide();
            this.inSide = new InSide();
            this.dip = Math.PI / 4;
            range = new double[4];
        }

        public Fornix(string name)
        {
            this.name = name;
            this.upFace = new UpFace();
            this.downFace = new DownFace();
            this.outSide = new OutSide();
            this.inSide = new InSide();
            this.dip = Math.PI / 4;
            range = new double[4];
        }

        /// <summary>
        /// 地层存入侧面上顶点
        /// </summary>
        /// <param name="vc"></param>
        public void createOutSideUpvers(VertexCollection vc)
        {
            this.outSide.addUpvers(vc);
        }
        /// <summary>
        /// 地层存入侧面下顶点
        /// </summary>
        /// <param name="vc"></param>
        public void createOutSideLowvers(VertexCollection vc)
        {
            this.outSide.addLowvers(vc);
        }
        /// <summary>
        /// 地层存入侧面下顶点
        /// </summary>
        public void createOutSideLowversFromSide()
        {
            for (int i = 0; i < this.outSide.sidelines.Count; ++i)
            {
                this.outSide.lowvers.addVer(this.outSide.sidelines[i].getVer(this.outSide.sidelines[i].Count - 1));
            }
        }

        /// <summary>
        /// 根据产状和贝塞尔系数生成几圈侧边点
        /// </summary>
        /// <param name="downElevation">底面海拔</param>
        /// <param name="Bessel_rad_cof">贝塞尔弧度系数</param>
        /// <param name="Bessel_num_control">贝塞尔控制点点数</param>
        /// <param name="Bessel_num_collected">贝塞尔采集点点数</param>
        public int createOutSideCircles(double downElevation, double Bessel_rad_cof, int Bessel_num_control, int Bessel_num_collected)
        {
            //由边界点及其产状生成侧边线控制点集合（边界点、中程控制点、底边点）
            createSidelinesControl(downElevation, Bessel_rad_cof, Bessel_num_control);
            //由侧边控制点生成Bessel_num_control + 1圈侧面顶点（边界点、贝塞尔曲线采集点、底边点），维持侧面顶点“位于同一侧边线”的关系
            return createSidelines(downElevation, Bessel_num_control, Bessel_num_collected);
        }
        
        /// <summary>
        /// 由边界点及其产状生成侧边线控制点集合（边界点、中程控制点、底边点）
        /// </summary>
        /// <param name="downElevation"></param>
        /// <param name="Bessel_rad_cof"></param>
        /// <param name="Bessel_num_control"></param>
        void createSidelinesControl(double downElevation, double Bessel_rad_cof, int Bessel_num_control)
        {
            Vertex vertexUpper, vertexLower;
            for (int i = 0; i < this.outSide.countUpVers(); ++i)
            {
                Sideline slc = new Sideline();//侧边控制点
                slc.addVer(this.outSide.getUpver(i));
                for (int j = 1; j < Bessel_num_control; ++j)
                {
                    vertexUpper = slc.getVer(j - 1);
                    vertexLower = new Vertex();
                    vertexLower.createControlVer(vertexUpper, downElevation, Bessel_rad_cof, Bessel_num_control);
                    slc.addVer(vertexLower);
                }
                this.outSide.sidelines_control.Add(slc);
            }
        }

        /// <summary>
        /// 由侧边控制点生成侧面顶点
        /// </summary>
        /// <param name="downElevation"></param>
        /// <param name="Bessel_num_control"></param>
        /// <param name="Bessel_num_collected"></param>
        /// <returns>总共新生成点的个数</returns>
        int createSidelines(double downElevation, int Bessel_num_control, int Bessel_num_collected)
        {
            Sideline curSLc, sl;
            for (int i = 0; i < this.outSide.sidelines_control.Count(); ++i)
            {
                curSLc = this.outSide.sidelines_control[i];
                sl = new Sideline();
                sl.createVerByBessel(curSLc, Bessel_num_collected, this.outSide.upvers.Count);
                this.outSide.sidelines.Add(sl);
            }
            return this.outSide.sidelines_control.Count() * (Bessel_num_collected - 1);
        }

        /// <summary>
        /// 地层贝塞尔侧面生成
        /// </summary>
        /// <returns></returns>
        public int createOutSideByBessel(int Bessel_num_collected)
        {
            int removeLowVersNum = 0;
            for (int i = 0; i < Bessel_num_collected - 1; ++i)
                removeLowVersNum += createOutSide(i);
            return removeLowVersNum;
        }

        /// <summary>
        /// 生成地层第i层侧面
        /// </summary>
        public int createOutSide(int sli)
        {
            //消除意外圈
            int removeLowVersNum = 0;
            removeLowVersNum += removeInsectedFace(sli);
            //消除狭小缝隙
            //removeLowVersNum += cancelMutateFace(sli);
            //侧面生成
            Triangle tri0, tri1;
            Vertex curV0, curV1,lowV0,lowV1;
            for (int i = 0; i < this.outSide.sidelines.Count; ++i)
            {
                curV0 = this.outSide.sidelines[i].getVer(sli);
                curV1 = i != this.outSide.sidelines.Count - 1?this.outSide.sidelines[i + 1].getVer(sli):this.outSide.sidelines[0].getVer(sli);
                lowV0 = this.outSide.sidelines[i].getVer(sli + 1);
                lowV1 = i != this.outSide.sidelines.Count - 1 ? this.outSide.sidelines[i + 1].getVer(sli + 1) : this.outSide.sidelines[0].getVer(sli + 1);
                tri0 = new Triangle();
                tri1 = new Triangle();
                tri0.addPoint(curV0);
                tri0.addPoint(curV1);
                tri0.addPoint(lowV0);
                this.outSide.addTri(tri0);
                if(!lowV1.toRemove)
                {
                    tri1.addPoint(curV1);
                    tri1.addPoint(lowV1);
                    tri1.addPoint(lowV0);
                    this.outSide.addTri(tri1);
                }
            }
            return removeLowVersNum;
        }

        /// <summary>
        /// 地层侧面生成
        /// </summary>
        public int createOutSide()
        {
            //消除意外圈
            int removeLowVersNum = 0;
            removeLowVersNum = removeInsectedFace();
            Triangle tri0, tri1;
            Vertex curV0, curV1;
            for (int i = 0; i < this.outSide.countUpVers(); ++i)
            {
                curV0 = this.outSide.getUpver(i);
                if (i != this.outSide.countUpVers()-1)
                    curV1 = this.outSide.getUpver(i + 1);
                else
                    curV1 = this.outSide.getUpver(0);
                if (curV0.corIdx == curV1.corIdx)//填入一个三角面
                {
                    tri0 = new Triangle();
                    tri0.addPoint(curV0);
                    tri0.addPoint(curV1);
                    tri0.addPoint(this.outSide.getDownver(curV0.corIdx));
                    this.outSide.addTri(tri0);
                }
                else//填入两个三角面
                {
                    tri0 = new Triangle();
                    tri1 = new Triangle();
                    tri0.addPoint(curV0);
                    tri0.addPoint(curV1);
                    tri0.addPoint(this.outSide.getDownver(curV0.corIdx));
                    tri1.addPoint(curV1);
                    tri1.addPoint(this.outSide.getDownver(curV1.corIdx));
                    tri1.addPoint(this.outSide.getDownver(curV0.corIdx));
                    this.outSide.addTri(tri0);
                    this.outSide.addTri(tri1);
                }
            }
            return removeLowVersNum;
        }

        /// <summary>
        /// 通过检查侧边底边线段有无相交（形成意外圈），来消除相交三角面，转而用其最外交点为顶点
        /// </summary>
        /// <returns>消除突变面导致底边边界点被移除的数量</returns>
        private int removeInsectedFace() 
        {
            int removeLowVersNum = 0;//移去侧面底边顶点个数
            int i, j, k = 0, idxEdge;//当前线段与其后第idxEdge个线段最后相交
            Vertex curV0, curV1, testV0, testV1, crossV0 = null, crossV1 = null;
            Vertex insertedVer;
            VertexCollection vcDown1 = new VertexCollection();

            for (i = 0; i < this.outSide.upvers.Count - 1; ++i)
            {
                curV0 = this.outSide.getDownver(i);
                curV1 = this.outSide.getDownver(i + 1);
                idxEdge = -1;
                //从当前线段后第二个线段开始，到最后一个线段，搜索交点（一个线段的下两个线段才可能与其“相交”，所以从i+2开始）
                for (j = i + 2; (i > 0 && j < this.outSide.upvers.Count) || (i == 0 && j < this.outSide.upvers.Count - 1); ++j)
                {
                    testV0 = this.outSide.getDownver(j);
                    if (j != this.outSide.lowvers.Count - 1)
                        testV1 = this.outSide.getDownver(j + 1);
                    else
                        testV1 = this.outSide.getDownver(0);
                    if (Tools.InsectionJudge(curV0, curV1, testV0, testV1))
                    {
                        crossV0 = testV0;
                        crossV1 = testV1;
                        idxEdge = j - i;
                    }
                }
                if (-1 == idxEdge)//说明没相交
                    continue;
                //“意外圈”中包含大多数边界点，意味着底边起始点可能处于“意外圈”中，
                //删除“意外圈”中所有顶点，起始点移至交点处
                //为防止这种情况，绘制地层边界时，就应当从地层界线外凸弧段部位开始绘制
                if (idxEdge > this.outSide.upvers.Count * 0.7)
                {
                    MessageBox.Show("底边起始点可能处于‘意外圈’中！");
                }
                removeLowVersNum += (idxEdge - 1);//移去idxEdge，又添加1个交点为新边界点
                //两个线段交点
                insertedVer = Tools.GetCrossPoint(curV0, curV1, crossV0, crossV1);
                insertedVer.ID = this.outSide.getDownver(i).ID + 1;
                //往新底面顶点集合中插入已确定的顶点
                for (j = k; j <= i; ++j)
                    vcDown1.addVer(this.outSide.lowvers.getVer(j));
                //更改上下顶点对应关系
                for (; j <= i + idxEdge; ++j)
                    this.outSide.getUpver(j).corIdx = vcDown1.Count - 1;
                //当前顶点后的所有顶点ID前移
                for (; j < this.outSide.lowvers.Count; ++j)
                {
                    this.outSide.getUpver(j).corIdx -= idxEdge;
                    this.outSide.getDownver(j).ID -= (idxEdge - 1);
                }
                vcDown1.addVer(insertedVer);
                k = i + idxEdge + 1;
                //向后推进到与当前侧面底边线段无交点的第一个线段
                i += idxEdge;
            }
            //插入剩余的顶点
            for (j = k; j < i + 1; ++j)
                vcDown1.addVer(this.outSide.lowvers.getVer(j));
            this.outSide.lowvers.clear();
            this.outSide.lowvers.addVerCollection(vcDown1);
            return removeLowVersNum;
        }

        /// <summary>
        /// 通过检查侧边第sli+1层线段有无相交（形成意外圈），来消除相交三角面，转而用其最外交点为顶点
        /// </summary>
        /// <param name="sli">侧面层数</param>
        /// <returns>消除突变面导致底边边界点被移除的数量</returns>
        private int removeInsectedFace(int sli)
        {
            int removeLowVersNum = 0;//移去侧面底边顶点个数
            int i, j, k = 0, idxEdge;//当前线段与其后第idxEdge个线段最后相交
            Vertex curV0, curV1, testV0, testV1, crossV0 = null, crossV1 = null;
            Vertex insertedVer;
            //VertexCollection vcDown1 = new VertexCollection();

            for (i = 0; i < this.outSide.upvers.Count - 1; ++i)
            {
                curV0 = this.outSide.sidelines[i].getVer(sli + 1);
                curV1 = this.outSide.sidelines[i+1].getVer(sli + 1);
                idxEdge = -1;
                //从当前线段后第二个线段开始，到最后一个线段，搜索交点（一个线段的下两个线段才可能与其“相交”，所以从i+2开始）
                for (j = i + 2; (i > 0 && j < this.outSide.sidelines.Count) || (i == 0 && j < this.outSide.sidelines.Count - 1); ++j)
                {
                    testV0 = this.outSide.sidelines[j].getVer(sli + 1);
                    if (j != this.outSide.sidelines.Count - 1)
                        testV1 = this.outSide.sidelines[j + 1].getVer(sli + 1);
                    else
                        testV1 = this.outSide.sidelines[0].getVer(sli + 1);
                    if (Tools.InsectionJudge(curV0, curV1, testV0, testV1))
                    {
                        crossV0 = testV0;
                        crossV1 = testV1;
                        idxEdge = j - i;
                    }
                }
                if (-1 == idxEdge)//说明没相交
                    continue;
                //“意外圈”中包含大多数边界点，意味着底边起始点可能处于“意外圈”中，
                //删除“意外圈”中所有顶点，起始点移至交点处
                //为防止这种情况，绘制地层边界时，就应当从地层界线外凸弧段部位开始绘制
                if (idxEdge > this.outSide.upvers.Count * 0.7)
                {
                    MessageBox.Show("底边起始点可能处于‘意外圈’中！");
                }
                //removeLowVersNum += (idxEdge - 1);//移去idxEdge，又添加1个交点为新边界点
                //两个线段交点
                insertedVer = Tools.GetCrossPoint(curV0, curV1, crossV0, crossV1);
                //insertedVer.ID = this.outSide.getDownver(i).ID + 1;
                //将i+1到i + idxEdge间的点移动到交点上
                for (j = i + 1; j <= i + idxEdge; ++j)
                {
                    this.outSide.sidelines[j].getVer(sli + 1).X(insertedVer.X());
                    this.outSide.sidelines[j].getVer(sli + 1).Y(insertedVer.Y());
                    this.outSide.sidelines[j].getVer(sli + 1).Z(insertedVer.Z());
                    if (j != i + 1)
                        this.outSide.sidelines[j].getVer(sli + 1).toRemove = true;
                }

                ////往新底面顶点集合中插入已确定的顶点
                //for (j = k; j <= i; ++j)
                //    vcDown1.addVer(this.outSide.lowvers.getVer(j));
                ////更改上下顶点对应关系
                //for (; j <= i + idxEdge; ++j)
                //    this.outSide.getUpver(j).corIdx = vcDown1.Count - 1;
                ////当前顶点后的所有顶点ID前移
                //for (; j < this.outSide.lowvers.Count; ++j)
                //{
                //    this.outSide.getUpver(j).corIdx -= idxEdge;
                //    this.outSide.getDownver(j).ID -= (idxEdge - 1);
                //}
                //vcDown1.addVer(insertedVer);
                k = i + idxEdge + 1;
                //向后推进到与当前侧面底边线段无交点的第一个线段
                i += idxEdge;
            }
            //插入剩余的顶点
            //for (j = k; j < i + 1; ++j)
            //    vcDown1.addVer(this.outSide.lowvers.getVer(j));
            //更新顶点集合
            //this.outSide.lowvers.clear();
            //this.outSide.lowvers.addVerCollection(vcDown1);

            return removeLowVersNum;
        }
    
        /// <summary>
        /// 解决突变面问题
        /// </summary>
        /// <param name="sli"></param>
        /// <returns></returns>
        private int cancelMutateFace(int sli)
        {
            int removeLowVersNum = 0;//移去侧面底边顶点个数
            int i, j, k = 0, idxVer;//当前线段与其后第idxEdge个线段最后相交
            Vertex curV0, curV1, testV, crossV0 = null, crossV1 = null;
            Vertex insertedVer;

            for (i = 0; i < this.outSide.upvers.Count - 1; ++i)
            {
                curV0 = this.outSide.sidelines[i].getVer(sli + 1);
                curV1 = this.outSide.sidelines[i + 1].getVer(sli + 1);
                testV = this.outSide.sidelines[i + 2].getVer(sli + 1);
                idxVer = 2;

                if (15.0 < Math.Acos(curV1.calAngle(curV0, testV))*180.0/Math.PI)
                    continue;

                //从当前线段后第二个线段开始，到最后一个线段，搜索交点（一个线段的下两个线段才可能与其“相交”，所以从i+2开始）
                for (j = i + 3; (i > 0 && j < this.outSide.sidelines.Count) || (i == 0 && j < this.outSide.sidelines.Count - 1); ++j)
                {
                    testV = this.outSide.sidelines[j].getVer(sli + 1);
                    ++idxVer;
                    if (15.0 < Math.Acos(curV1.calAngle(curV0, testV)) * 180.0 / Math.PI)
                        break;
                }
                //若缝合线长度大于curV01的1.5倍则不做处理
                if (curV1.calDistance(testV) > curV1.calDistance(curV0) * 1.5)
                    continue;

                if (-1 == idxVer)//说明没相交
                    continue;
                //“意外圈”中包含大多数边界点，意味着底边起始点可能处于“意外圈”中，
                //删除“意外圈”中所有顶点，起始点移至交点处
                //为防止这种情况，绘制地层边界时，就应当从地层界线外凸弧段部位开始绘制
                if (idxVer > this.outSide.upvers.Count * 0.7)
                {
                    MessageBox.Show("底边起始点可能处于‘意外圈’中！");
                }
                //removeLowVersNum += (idxEdge - 1);//移去idxEdge，又添加1个交点为新边界点
                //两个线段交点
                insertedVer = Tools.GetCrossPoint(curV0, curV1, crossV0, crossV1);
                //insertedVer.ID = this.outSide.getDownver(i).ID + 1;
                //将i+1到i + idxEdge间的点移动到交点上
                for (j = i + 2; j <= i + idxVer; ++j)
                {
                    this.outSide.sidelines[j].getVer(sli + 1).X(insertedVer.X());
                    this.outSide.sidelines[j].getVer(sli + 1).Y(insertedVer.Y());
                    this.outSide.sidelines[j].getVer(sli + 1).Z(insertedVer.Z());
                    if (j != i + 1)
                        this.outSide.sidelines[j].getVer(sli + 1).toRemove = true;
                }
                k = i + idxVer + 1;
                //向后推进到与当前侧面底边线段无交点的第一个线段
                i += idxVer;
            }
            return removeLowVersNum;
        }

        public void createInSide(Fornix foxnix)
        {
            this.inSide.createInByOut(foxnix);
        }

        public void createUpFace(VertexCollection rasterPoints)
        {
            this.upFace.createBySide(this.outSide.upvers, this.inSide.upvers, rasterPoints);
        }

        public void createDownFace()
        {
            this.downFace.createBySide(this.outSide.lowvers, this.inSide.lowvers);
        }

        /// <summary>
        /// 设置地层范围
        /// </summary>
        /// <param name="xmin"></param>
        /// <param name="ymin"></param>
        /// <param name="xmax"></param>
        /// <param name="ymax"></param>
        public void setRange(double xmin, double ymin, double xmax, double ymax)
        {
            range[0] = xmin;
            range[1] = ymin;
            range[2] = xmax;
            range[3] = ymax;
        }

        /// <summary>
        /// 对比地层和等高线范围
        /// </summary>
        /// <param name="contour">等高线</param>
        /// <returns>是否可以快速排除</returns>
        public bool fastJudge(Contour contour)
        {
            if(this.range[0] > contour.range[2]
                || this.range[1] > contour.range[3]
                || this.range[2] < contour.range[0]
                || this.range[3] < contour.range[1])
                return true;
            return false;
        }

    }

    class UpFace
    {
        public SortedList<string, Triangle> tris { get; private set; }
        public VertexCollection InnerPoints { get; private set; }

        public UpFace()
        {
            this.tris = new SortedList<string, Triangle>();
        }
        public void createBySide(VertexCollection outSideLowvers, VertexCollection inSideLowvers, VertexCollection pointsOfDEM)
        {
            this.InnerPoints = Delaunay.GenTriMesh(tris, outSideLowvers, inSideLowvers, pointsOfDEM);
        }
    }

    class DownFace
    {
        public SortedList<string, Triangle> tris { get; private set; }

        public DownFace()
        {
            this.tris = new SortedList<string, Triangle>();
        }

        public void createBySide(VertexCollection outSideLowvers, VertexCollection inSideLowvers)
        {

            Delaunay.GenTriMesh(tris, outSideLowvers, inSideLowvers, null);
        }
    }

    class OutSide
    {
        public SortedList<string, Triangle> tris { get; private set; }
        public List<Edge> edges { get; private set; }
        public VertexCollection upvers { get; private set; }
        public VertexCollection lowvers { get; private set; }
        public List<Sideline> sidelines_control;//侧边线控制点集合
        public List<Sideline> sidelines;//侧边线顶点集合
        public OutSide()
        {
            this.tris = new SortedList<string, Triangle>();
            this.edges = new List<Edge>();
            this.upvers = new VertexCollection();
            this.lowvers = new VertexCollection();
            this.sidelines_control = new List<Sideline>();
            this.sidelines = new List<Sideline>();
        }
        public void addUpvers(VertexCollection upvers)
        {
            this.upvers.addVerCollection(upvers);
        }
        public void addLowvers(VertexCollection lowvers)
        {
            this.lowvers.addVerCollection(lowvers);
        }
        public Vertex getUpver(int index)
        {
            return this.upvers.getVer(index);
        }
        public Vertex getDownver(int index)
        {
            return this.lowvers.getVer(index);
        }
        public int countUpVers()
        {
            return this.upvers.Count;
        }
        public void addTri(Triangle tri)
        {
            this.tris.Add(tri.name, tri);
        }
        public void addEdge(Edge edge)
        {
            this.edges.Add(edge);
        }
    }
    class InSide
    {
        private SortedList<string, Triangle> tris;
        public VertexCollection upvers { get; private set; }
        public VertexCollection lowvers { get; private set; }

        public InSide()
        {
            this.tris = new SortedList<string, Triangle>();
            this.upvers = new VertexCollection();
            this.lowvers = new VertexCollection();
        }

        public void createInByOut(Fornix fornix)
        {
            this.upvers.clear();
            this.lowvers.clear();
            this.upvers.addVerCollection(fornix.outSide.upvers);
            this.lowvers.addVerCollection(fornix.outSide.lowvers);
        }
    }

    class Triangle
    {
        public string name { get; set; }
        public VertexCollection points { get; private set; }  //三角面点集
        public bool isCavity { get; set; }

        public Triangle()
        {
            this.points = new VertexCollection();
            this.points.clear();
            this.isCavity = false;
        }

        public void addPoint(Vertex ver)
        {
            this.points.addVer(ver);
            if (this.points.Count == 3)
            {
                string point0 = this.points.getVer(0).ID.ToString(),
                    point1 = this.points.getVer(1).ID.ToString(),
                    point2 = this.points.getVer(2).ID.ToString();
                if (this.points.getVer(0).innerPoint)
                    point0 = "in" + point0;
                if (this.points.getVer(1).innerPoint)
                    point1 = "in" + point1;
                if (this.points.getVer(2).innerPoint)
                    point2 = "in" + point2;
                this.name = point0 + "_" + point1 + "_" + point2;
            }
        }

        public Vertex calCircumCenter()
        {
            if (this.points.Count != 3)
                MessageBox.Show("顶点不为3，无法求外心！");
            return this.points.getVer(0).calCircumCenter(this.points.getVer(1), this.points.getVer(2));
        }

        public List<Triangle> getAdjTri(SortedList<string, Triangle> tris)
        {
            Vertex ver0 = this.points.getVer(0);
            Vertex ver1 = this.points.getVer(1);
            Vertex ver2 = this.points.getVer(2);
            List<Triangle> adjTris = new List<Triangle>();
            adjTris.Clear();
            for (int i = 0; i < tris.Count; ++i)
            {
                Triangle tri = tris.Values[i];
                if ((tri.hasVer(ver0) && tri.hasVer(ver1) && !tri.hasVer(ver2))
                    || (tri.hasVer(ver0) && !tri.hasVer(ver1) && tri.hasVer(ver2))
                    || (!tri.hasVer(ver0) && tri.hasVer(ver1) && tri.hasVer(ver2)))
                    adjTris.Add(tri);
            }
            return adjTris;
        }

        public bool hasVer(Vertex ver)
        {
            for (int i = 0; i < 3; ++i)
                if (ver.ID == this.points.getVer(i).ID && !(ver.innerPoint ^ this.points.getVer(i).innerPoint))
                    return true;
            return false;
        }

        public bool hasVer(int ID)
        {
            for (int i = 0; i < 3; ++i)
                if (ID == this.points.getVer(i).ID)
                    return true;
            return false;
        }

        /// <summary>
        /// 判断每个顶点都在环上
        /// </summary>
        /// <param name="minID"></param>
        /// <param name="maxID"></param>
        /// <returns></returns>
        public bool pointsOnRing(int minID, int maxID)
        {
            if ((this.points.getVer(0).ID >= minID && this.points.getVer(0).ID <= maxID)
                && (this.points.getVer(1).ID >= minID && this.points.getVer(1).ID <= maxID)
                && (this.points.getVer(2).ID >= minID && this.points.getVer(2).ID <= maxID))
                return true;
            return false;
        }

        public double calVector()
        {
            return this.points.getVer(1).calVector(this.points.getVer(0), this.points.getVer(2));
        }

        /// <summary>
        /// 按照三角形中点ID从小到大求
        /// </summary>
        /// <param name="featureLayer"></param>
        /// <param name="fornixs"></param>
        /// <returns></returns>
        public double calOrderVector()
        {
            return this.points.getVer(1).calVector(this.points.getVer(0), this.points.getVer(2));
        }

        public Triangle reverse()
        {
            List<Vertex> revs = new List<Vertex>();

            revs.Add(this.points.getVer(2));
            revs.Add(this.points.getVer(1));
            revs.Add(this.points.getVer(0));

            this.points.remove(2);
            this.points.remove(1);
            this.points.remove(0);

            this.addPoint(revs[0]);
            this.addPoint(revs[1]);
            this.addPoint(revs[2]);

            return this;
        }

        public string reverseName()
        {
            string point0 = this.points.getVer(0).ID.ToString(),
                   point1 = this.points.getVer(1).ID.ToString(),
                   point2 = this.points.getVer(2).ID.ToString();
            if (this.points.getVer(0).innerPoint)
                point0 = "in" + point0;
            if (this.points.getVer(1).innerPoint)
                point1 = "in" + point1;
            if (this.points.getVer(2).innerPoint)
                point2 = "in" + point2;
            return point2 + "_" + point1 + "_" + point0;
        }

    }

    class Edge
    {
        public Vertex v0 { get; set; }
        public Vertex v1 { get; set; }
        public String name;
        public String[] adjTriName { get; set; }//删除底边相交点时，用于找到应当删除的三角面

        public Edge(Vertex v0, Vertex v1)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.adjTriName = new string[2];
            name = v0.ID.ToString() + "_" + v1.ID.ToString();
        }

        public bool Equals(Edge edge)
        {
            if (this.v0.ID == edge.v0.ID && this.v1.ID == edge.v1.ID
                || this.v0.ID == edge.v1.ID && this.v1.ID == edge.v0.ID)
                return true;
            return false;
        }

        public bool Equals(Vertex v0, Vertex v1)
        {
            if ((this.v0.ID == v0.ID && this.v1.ID == v1.ID && !(this.v0.innerPoint ^ v0.innerPoint) && !(this.v1.innerPoint ^ v1.innerPoint))
                || (this.v0.ID == v1.ID && this.v1.ID == v0.ID && !(this.v0.innerPoint ^ v1.innerPoint) && !(this.v1.innerPoint ^ v0.innerPoint)))
                return true;
            return false;
        }

        public Vertex getVer(int idx)
        {
            if (0 == idx)
                return this.v0;
            else
                return this.v1;
        }

        /// <summary>
        /// 判断俩直线是否相交
        /// </summary>
        /// <param name="edge">第二条线段</param>
        /// <returns>是否相交</returns>
        public bool InsectionJudge(Edge edge)
        {
            Vertex a = this.v0;
            Vertex b = this.v1;
            Vertex c = edge.v0;
            Vertex d = edge.v1;
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
        public Vertex getCrossPoint(Edge edge)
        {
            if (null == edge)
                return null;
            Vertex crossPoint = new Vertex();
            Vertex a = this.v0;
            Vertex b = this.v1;
            Vertex c = edge.v0;
            Vertex d = edge.v1;
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

    class Contour :
        VertexCollection
    {
        public double[] range { get; private set; }

        public Contour()
        {
            range = new double[4];
        }
        public void setRange(double xmin, double ymin, double xmax, double ymax)
        {
            range[0] = xmin;
            range[1] = ymin;
            range[2] = xmax;
            range[3] = ymax;
        }
    }

    class VertexCollection
    {
        private List<Vertex> vers;
        public int Count { get { return this.vers.Count; } }

        public VertexCollection()
        {
            this.vers = new List<Vertex>();
            this.clear();
        }
        public void addVer(Vertex ver)
        {
            this.vers.Add(ver);
        }
        public void addVer(Vertex ver, int idx)
        {
            this.vers.Insert(idx, ver);
        }
        public void addVerCollection(VertexCollection vers)
        {
            int num = vers.Count;
            for (int i = 0; i < num; ++i)
            {
                this.addVer(vers.getVer(i));
            }
        }
        public void clear()
        {
            this.vers.Clear();
        }
        public Vertex getVer(int index)
        {
            return this.vers[index];
        }
        public void remove(int index)
        {
            this.vers.Remove(this.vers[index]);
        }

        public void removeVers(int index, int count)
        {
            this.vers.RemoveRange(index, count);
        }

        public void insert(int idx, Vertex ver)
        {
            this.vers.Insert(idx, ver);
        }
    }

    /// <summary>
    /// 侧边线顶点集合/侧边线控制点
    /// </summary>
    class Sideline : VertexCollection
    {
        /// <summary>
        /// 由控制点生成贝塞尔曲线，得到采集点
        /// </summary>
        /// <param name="slc"></param>
        /// <param name="num"></param>
        /// <param name="num_upvers">边界点个数</param>
        public void createVerByBessel(Sideline slc,int num,int num_upvers)
        {
            Vertex verBessel;
            this.addVer(slc.getVer(0));//插入顶边点
            double t = 0;
            t += (1.0 / (num - 1.0));
            int i;
            for (i = 1; i < num - 1; ++i)
            {
                verBessel = calBesselVer(slc, t);
                if (null != verBessel)
                {
                    verBessel.ID = slc.getVer(0).ID + num_upvers * i;
                    this.addVer(verBessel);
                }
                t += (1.0 / (num - 1.0));
            }
            slc.getVer(slc.Count - 1).ID = slc.getVer(0).ID + num_upvers * i;
            this.addVer(slc.getVer(slc.Count - 1));//插入底边点
        }

        /// <summary>
        /// 计算贝塞尔曲线采集点坐标
        /// </summary>
        /// <param name="slc"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        private Vertex calBesselVer(Sideline slc, double t)
        {
            if (3 == slc.Count)//一个中程控制点
            {
                Vertex verBessel = new Vertex();
                verBessel.X(Math.Pow(1.0 - t, 2) * slc.getVer(0).X() + 2 * t * (1 - t) * slc.getVer(1).X() + t * t * slc.getVer(2).X());
                verBessel.Y(Math.Pow(1.0 - t, 2) * slc.getVer(0).Y() + 2 * t * (1 - t) * slc.getVer(1).Y() + t * t * slc.getVer(2).Y());
                verBessel.Z(Math.Pow(1.0 - t, 2) * slc.getVer(0).Z() + 2 * t * (1 - t) * slc.getVer(1).Z() + t * t * slc.getVer(2).Z());
                return verBessel;
            }
            else if (4 == slc.Count)//两个中程控制点
            {
                //之后再写
            }
            else if (4 < slc.Count)
                MessageBox.Show("暂不考虑高阶贝塞尔曲线");
            return null;
        }
    }


    /// <summary>
    /// 地层等高线交点
    /// </summary>
    class InsectVer : Vertex
    {
        public int VerIdx { get; set; }//上邻接边界点ID
    }

    class Vertex
    {
        public int ID { get; set; }
        public int corIdx { get; set; }//底边/顶边对应边界点在其集合内的编号（不是ID）
        public List<int> sideedgeVerIDs;//侧边线内部顶点（贝塞尔）
        private Coordinate position { get; set; }
        public Occurrence occurrence { get; set; }
        public bool innerPoint { get; set; }//是否是内部点
        public bool toRemove { get; set; }//消除自相交或狭小缝隙的时候有没有被移除
        public int type { get; set; }//类型：0-边界点 1-底边边界点 2-等高线折点
        public double angle{ get; set; }//角度
        public double curvature { get; set; }//曲率
        public Vertex()
        {
            this.position = new Coordinate();
            this.occurrence = new Occurrence();
            this.innerPoint = false;
            this.toRemove = false;
        }

        /// <summary>
        /// 计算顶点产状(外接圆法)
        /// </summary>
        /// <param name="prev"></param>
        /// <param name="nextv"></param>
        public void calOccuurence(Vertex prev, Vertex nextv)
        {
            Vertex circumCenter = calCircumCenter(prev, nextv);
            //外心指向当前点的向量
            double x = this.X() - circumCenter.X();
            double y = this.Y() - circumCenter.Y();
            if (x > 0 && y > 0)
            {
                this.occurrence.inclination = Math.Atan(Math.Abs(x) / Math.Abs(y));
            }
            else if (x > 0 && y < 0)
            {
                this.occurrence.inclination = Math.PI / 2 + Math.Atan(Math.Abs(y) / Math.Abs(x));
            }
            else if (x < 0 && y < 0)
            {
                this.occurrence.inclination = Math.PI + Math.Atan(Math.Abs(x) / Math.Abs(y));
            }
            else
            {
                this.occurrence.inclination = Math.PI * 3 / 2 + Math.Atan(Math.Abs(y) / Math.Abs(x));
            }
            //规定倾角指向地层线外
            if (0.0 < this.calVector(prev, nextv))
            {
                this.occurrence.inclination = (this.occurrence.inclination - Math.PI >= 0.0) ?
                    this.occurrence.inclination - Math.PI : this.occurrence.inclination + Math.PI;
            }
            //计算走向
            this.occurrence.trend = (this.occurrence.inclination - Math.PI / 2 >= 0.0) ?
                this.occurrence.inclination - Math.PI / 2 : this.occurrence.inclination + Math.PI / 2;
        }

        /// <summary>
        /// 计算顶点产状(法向量法)
        /// </summary>
        /// <param name="prev"></param>
        /// <param name="nextv"></param>
        public void calOccuurence1(Vertex prev, Vertex nextv)
        {
            double A = (this.Y() - prev.Y()) * (nextv.Z() - prev.Z()) - (this.Z() - prev.Z()) * (nextv.Y() - prev.Y());
            double B = (this.Z() - prev.Z()) * (nextv.X() - prev.X()) - (this.X() - prev.X()) * (nextv.Z() - prev.Z());
            double C = (this.X() - prev.X()) * (nextv.Y() - prev.Y()) - (this.Y() - prev.Y()) * (nextv.X() - prev.X());
            double row1 = Math.Atan(Math.Abs(B) / Math.Abs(A));
            if (A > 0 && B > 0)
                this.occurrence.inclination = Math.PI / 2 - row1;
            else if (A < 0 && B > 0)
                this.occurrence.inclination = Math.PI * 3 / 2 + row1;
            else if (A < 0 && B < 0)
                this.occurrence.inclination = Math.PI * 3 / 2 - row1;
            else if (A > 0 && B < 0)
                this.occurrence.inclination = Math.PI / 2 + row1;

            //规定倾角指向地层线外
            double Xpt = this.X() - prev.X();
            double Ypt = this.Y() - prev.Y();
            double Convexity = this.calVector(prev, nextv);
            double Angle = Math.Atan(this.calAngle(prev, nextv));
            double Angle_prev_normal = Math.Atan(((-Xpt) * A + (-Ypt) * B) / (Math.Sqrt(Xpt * Xpt + Ypt * Ypt) + Math.Sqrt(A * A + B * B)));
            //计算倾向
            if (0.0 > Convexity)
                if (0.0 > Xpt * B - Ypt * A && Angle > Angle_prev_normal)
                    this.occurrence.inclination = (this.occurrence.inclination - Math.PI >= 0.0) ?
                    this.occurrence.inclination - Math.PI : this.occurrence.inclination + Math.PI;
            else
            {
                if (0.0 > Xpt * B - Ypt * A)
                {
                    if (Math.PI > Angle_prev_normal+Angle)
                        //外接圆法
                        this.calOccuurence(prev, nextv);
                    else
                        this.occurrence.inclination = (this.occurrence.inclination - Math.PI >= 0.0) ?
                    this.occurrence.inclination - Math.PI : this.occurrence.inclination + Math.PI;
                }
                else if (Angle_prev_normal > Angle) 
                    //外接圆法
                    this.calOccuurence(prev, nextv);
            }

            //计算走向
            this.occurrence.trend = (this.occurrence.inclination - Math.PI / 2 >= 0.0) ?
                this.occurrence.inclination - Math.PI / 2 : this.occurrence.inclination + Math.PI / 2;

            //三点法计算倾角
            this.occurrence.dip = Math.Acos(Math.Abs(C) / Math.Sqrt(A * A + B * B + C * C));
            if (double.IsNaN(this.occurrence.dip))
                Console.WriteLine("");
        }

        /// <summary>
        /// 根据相邻点修复偏得比较厉害的倾角
        /// </summary>
        public void fixInclination()
        {
            
        }

        public double calDip(Vertex prev, Vertex nextv)
        {
            double A = (this.Y() - prev.Y()) * (nextv.Z() - prev.Z()) - (this.Z() - prev.Z()) * (nextv.Y() - prev.Y());
            double B = (this.Z() - prev.Z()) * (nextv.X() - prev.X()) - (this.X() - prev.X()) * (nextv.Z() - prev.Z());
            double C = (this.X() - prev.X()) * (nextv.Y() - prev.Y()) - (this.Y() - prev.Y()) * (nextv.X() - prev.X());
            return Math.Acos(Math.Abs(C) / Math.Sqrt(Math.Pow(A, 2) + Math.Pow(B, 2) + Math.Pow(C, 2)));
        }

        public Vertex calCircumCenter(Vertex prev, Vertex nextv)
        {
            Vertex circumCenter = new Vertex();
            double A0 = prev.X() * prev.X() + prev.Y() * prev.Y();
            double A1 = this.X() * this.X() + this.Y() * this.Y();
            double A2 = nextv.X() * nextv.X() + nextv.Y() * nextv.Y();
            double Ax = A1 * nextv.Y() + A0 * this.Y() + A2 * prev.Y() - prev.Y() * A1 - A0 * nextv.Y() - this.Y() * A2;
            double Ay = this.X() * A2 + prev.X() * A1 + A0 * nextv.X() - A0 * this.X() - prev.X() * A2 - A1 * nextv.X();
            double B = this.X() * nextv.Y() + prev.X() * this.Y() + prev.Y() * nextv.X() - prev.Y() * this.X() - prev.X() * nextv.Y() - this.Y() * nextv.X();
            if (B == 0.0)
                MessageBox.Show("边界点可能重复！");
            circumCenter.X(Ax / (2 * B));
            circumCenter.Y(Ay / (2 * B));
            return circumCenter;
        }

        /// <summary>
        /// 创建顶面顶点的底面对应顶点（厚度）
        /// </summary>
        /// <param name="ver"></param>
        /// <returns></returns>
        public Vertex createDownVer(Vertex ver)
        {
            double height = 100.0;
            this.X(ver.X() + height * Math.Sin(ver.occurrence.inclination) / Math.Tan(ver.occurrence.dip));
            this.Y(ver.Y() + height * Math.Cos(ver.occurrence.inclination) / Math.Tan(ver.occurrence.dip));
            this.Z(ver.Z() - height);
            return this;
        }

        /// <summary>
        /// 创建顶面顶点的底面对应顶点（底面海拔）
        /// </summary>
        /// <param name="ver"></param>
        /// <param name="downElevation"></param>
        /// <returns></returns>
        public Vertex createDownVer(Vertex ver, double downElevation)
        {
            double height = ver.Z() - downElevation;
            this.X(ver.X() + height * Math.Sin(ver.occurrence.inclination) / Math.Tan(ver.occurrence.dip));
            this.Y(ver.Y() + height * Math.Cos(ver.occurrence.inclination) / Math.Tan(ver.occurrence.dip));
            this.Z(downElevation);
            if (double.IsNaN(this.X()) || double.IsNaN(this.Y()) || double.IsNaN(this.Z()))
                Console.WriteLine("");
            return this;
        }

        /// <summary>
        /// 由边界点计算侧边控制点
        /// </summary>
        /// <param name="ver"></param>
        /// <param name="downElevation"></param>
        /// <param name="Bessel_rad_cof"></param>
        /// <param name="Bessel_num_control"></param>
        /// <returns></returns>
        public Vertex createControlVer(Vertex ver, double downElevation, double Bessel_rad_cof, int Bessel_num_control)
        {
            double height = (ver.Z() - downElevation) / (Bessel_num_control - 1);
            this.X(ver.X() + height * Math.Sin(ver.occurrence.inclination) / Math.Tan(ver.occurrence.dip));
            this.Y(ver.Y() + height * Math.Cos(ver.occurrence.inclination) / Math.Tan(ver.occurrence.dip));
            this.Z(ver.Z() + (downElevation - ver.Z()) / (Bessel_num_control - 1));
            this.occurrence.inclination = ver.occurrence.inclination;
            this.occurrence.dip = ver.occurrence.dip * Bessel_rad_cof;
            if (double.IsNaN(this.X()) || double.IsNaN(this.Y()) || double.IsNaN(this.Z()))
                Console.WriteLine("");
            return this;
        }

        /// <summary>
        /// 计算夹角
        /// </summary>
        /// <param name="prev"></param>
        /// <param name="nextv"></param>
        /// <returns></returns>
        public double calAngle(Vertex prev, Vertex nextv)
        {
            double x0 = prev.X() - this.X();
            double y0 = prev.Y() - this.Y();
            double x1 = nextv.X() - this.X();
            double y1 = nextv.Y() - this.Y();
            return (x0 * x1 + y0 * y1) / (Math.Sqrt(x0 * x0 + y0 * y0) * Math.Sqrt(x1 * x1 + y1 * y1));
        }

        /// <summary>
        /// 计算外积
        /// </summary>
        /// <param name="prev"></param>
        /// <param name="nextv"></param>
        /// <returns></returns>
        public double calVector(Vertex prev, Vertex nextv)
        {
            double x0 = this.X() - prev.X();
            double y0 = this.Y() - prev.Y();
            double x1 = nextv.X() - this.X();
            double y1 = nextv.Y() - this.Y();
            return x0 * y1 - x1 * y0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="prev"></param>
        /// <param name="nextv"></param>
        /// <returns></returns>
        public double calCurvature(Vertex prev, Vertex nextv)
        {
            this.curvature = this.angle / (this.calDistance(prev) + this.calDistance(nextv));
            return this.curvature;
        }

        public double calDistance(Vertex ver)
        {
            return Math.Sqrt(Math.Pow(this.X() - ver.X(), 2) + Math.Pow(this.Y() - ver.Y(), 2));
        }

        public bool inside(Triangle tri)
        {
            double CAd = tri.points.getVer(0).calVector(tri.points.getVer(2), this);
            double dAB = tri.points.getVer(0).calVector(this, tri.points.getVer(1));
            double ABd = tri.points.getVer(1).calVector(tri.points.getVer(0), this);
            double dBC = tri.points.getVer(1).calVector(this, tri.points.getVer(2));
            if (CAd * dAB > 0.0 && ABd * dBC > 0.0)
                return true;
            return false;
        }

        public double X() { return this.position.x; }
        public double Y() { return this.position.y; }
        public double Z() { return this.position.z; }
        public void X(double x) { this.position.x = x; }
        public void Y(double y) { this.position.y = y; }
        public void Z(double z) { this.position.z = z; }

    }

    class Occurrence
    {
        public double trend { get; set; }//走向
        public double inclination { get; set; }//倾向
        public double dip { get; set; }//倾角
    }

    class Coordinate
    {
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
    }

    class Tools
    {
        /// <summary>
        /// 判断俩直线是否相交
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        /// <returns>是否相交</returns>
        public static bool InsectionJudge(Vertex a,Vertex b,Vertex c,Vertex d)
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
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        /// <returns>交点</returns>
        public static InsectVer GetCrossPoint(Vertex a, Vertex b, Vertex c, Vertex d)
        {
            InsectVer crossPoint = new InsectVer();
            if (b.Y() != a.Y() && d.Y() != c.Y())
            {
                double A = (b.X() - a.X()) / (b.Y() - a.Y());
                double B = (d.X() - c.X()) / (d.Y() - c.Y());
                double y = (c.X() - a.X() + A * a.Y() - B * c.Y()) / (A - B);
                double x = A * y - A * a.Y() + a.X();
                crossPoint.X(x);
                crossPoint.Y(y);
                if (double.IsNaN(x) || double.IsNaN(y))
                    Console.Write("");
                crossPoint.Z((a.Z() + b.Z() + c.Z() + d.Z()) / 4);
                crossPoint.ID = b.ID;
                return crossPoint;
            }
            else if (b.Y() == a.Y())
            {
                crossPoint.X(c.X() + (d.X() - c.X()) * (a.Y() - c.Y()) / (d.Y() - c.Y()));
                crossPoint.Y(a.Y());
                return crossPoint;
            }
            else if (d.Y() == c.Y())
            {
                crossPoint.X(c.X() + (b.X() - a.X()) * (c.Y() - a.Y()) / (b.Y() - a.Y()));
                crossPoint.Y(c.Y());
                return crossPoint;
            }
            else
                return null;
            
        }

    }
}
