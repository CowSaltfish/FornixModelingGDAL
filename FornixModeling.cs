using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

using OSGeo.OGR;
//using OSGeo.OSR;
//using OSGeo.GDAL;

namespace FornixModelingGDAL
{
    class FornixModeling
    {
        private static string outpath;
        private static string imgpath;
        private static string shppath;
        private static DataSource _ds;
        private static List<Fornix> _fornixs = new List<Fornix>();//穹窿地层集合
        private static int _verNum = 1;
        private static List<Layer> _layers = new List<Layer>();
        private static Layer _contourLayer;
        private static List<Contour> _contours = new List<Contour>();
        private static List<VertexCollection> _Intersects = new List<VertexCollection>();//地层界线与等高线交点集合
        private static VertexCollection _RasterPoints;
        private static double _downElevation;//当按照底面海拔建模时，底面海拔
        private static double Bessel_rad_cof;//贝塞尔弧度系数
        private static int Bessel_num_control;//贝塞尔控制点数
        private static int Bessel_num_collected;//贝塞尔采集点数

        /// <summary>
        /// 生成穹窿模型
        /// </summary>
        public static void GenModel()
        {
            //Gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
            //Gdal.SetConfigOption("SHAPE_ENCODING", "");

            outpath = @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\Huangling\Export\HLQL1.obj";
            imgpath = @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\Huangling\DEM\input_dem.img";
            shppath = @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\Huangling\Export\shp\";
            _downElevation = -300.0;
            Bessel_num_control = 3;
            Bessel_num_collected = 4;
            Bessel_rad_cof = 1.5;
            //由带值点要素生成穹窿模型
            createFornixModelByPoint();
            MessageBox.Show("穹窿生成成功！");
        }

        /// <summary>
        /// 由带值点要素生成穹窿模型
        /// </summary>
        private static void createFornixModelByPoint()
        {
            //采样点
            //_contour = DEMServer.CalContourByGdal(vctPath, 300);
            _RasterPoints = DEMServer.getVersFromDEM(imgpath);

            Fornix preFornix = null;

            //灵岩山矢量数据
            //ReadShp(0, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\LingYanShan\codedata\vertex1\lyspv1.shp");
            //ReadShp(0, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\LingYanShan\codedata\vertex1\lyspv2.shp");
            //ReadShp(0, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\LingYanShan\codedata\vertex1\lyspv3.shp");

            //黄陵穹窿矢量数据
            ReadShp(0, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\Huangling\Vector\Strata11_Points.shp");
            ReadShp(0, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\Huangling\Vector\Strata10_Points.shp");
            //等高线折点
            ReadShp(1, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\Huangling\Vector\contour\FeatureVerticesToPoi\contour300_points.shp");
            //产状点(先模拟)
            //ReadShp(2, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\LingYanShan\codedata\vertex1\产状点.shp");

            //将等高线折点存入集合
            createContour();
            int si = 0;
            foreach (Layer PointLayer in _layers)
            {
                Fornix fornix = new Fornix(PointLayer.GetName());
                fornix.Bessel_num_collected = Bessel_num_collected;
                //读入一个地层
                createStrata(PointLayer, fornix, si);
                //计算地层边界点角度与曲率
                //calVertexAC(fornix);
                //计算地层边界点产状
                calOccurrence(fornix);
                //将顶面边界点产状信息插回shp中
                //SetShp(PointLayer, fornix);
                //侧面生成
                createSide(fornix, preFornix);
                //组合
                _fornixs.Add(fornix);
                preFornix = fornix;
                ++si;
            }
            foreach (Fornix fornix in _fornixs)
            {
                //基于相邻顶边，生成顶面
                fornix.createUpFace(_RasterPoints);
                //基于相邻底边，生成底面
                fornix.createDownFace();
            }
            //输出模型底面shp
            writeDownSHP(_fornixs, shppath);
            //打印obj
            ObjWriter.writeFornixObj(_fornixs, outpath, _verNum);
        }

        /// <summary>
        /// 将等高线折点存入集合
        /// </summary>
        private static void createContour()
        {
            double curContour = _contourLayer.GetFeature(0).GetFieldAsDouble("CONTOUR");
            Contour contourPoints = new Contour();
            double xmin = 10e9, ymin = 10e9, xmax = -10e9, ymax = -10e9;
            for (int i = 0; i < (int)_contourLayer.GetFeatureCount(0); ++i)
            {
                Vertex contourPoint = new Vertex();
                contourPoint.X(_contourLayer.GetFeature(i).GetFieldAsDouble("POINT_X"));
                if (xmin > contourPoint.X())
                    xmin = contourPoint.X();
                if (xmax < contourPoint.X())
                    xmax = contourPoint.X();
                contourPoint.Y(_contourLayer.GetFeature(i).GetFieldAsDouble("POINT_Y"));
                if (ymin > contourPoint.Y())
                    ymin = contourPoint.Y();
                if (ymax < contourPoint.Y())
                    ymax = contourPoint.Y();
                contourPoint.Z(_contourLayer.GetFeature(i).GetFieldAsDouble("CONTOUR"));
                contourPoint.type = 2;
                if (curContour == contourPoint.Z())
                    contourPoints.addVer(contourPoint);
                else
                {
                    contourPoints.setRange(xmin, ymin, xmax, ymax);
                    _contours.Add(contourPoints);
                    contourPoints = new Contour();
                    curContour = contourPoint.Z();
                    contourPoints.addVer(contourPoint);
                    xmin = 10e9; ymin = 10e9; xmax = -10e9; ymax = -10e9;
                }
            }
        }

        /// <summary>
        /// 读入一个地层
        /// </summary>
        /// <param name="PointLayer">地层边界点图层</param>
        /// <param name="fornix">地层容器</param>
        /// <param name="si">第几个地层（仅仅是为了按照一定的规律赋予几个边界点随机的倾角，模拟）</param>
        private static void createStrata(Layer PointLayer, Fornix fornix,int si)
        {
            VertexCollection vcUp = fornix.outSide.upvers;
            Random rd = new Random();
            double xmin = 10e9, ymin = 10e9, xmax = -10e9, ymax = -10e9;
            for (int i = 0; i < (int)PointLayer.GetFeatureCount(0) - 1; ++i)
            {
                Vertex vertex = new Vertex();
                vertex.ID = _verNum++;
                //PointLayer.GetFeature(i).
                vertex.X(PointLayer.GetFeature(i).GetFieldAsDouble("POINT_X"));
                if (xmin > vertex.X())
                    xmin = vertex.X();
                if (xmax < vertex.X())
                    xmax = vertex.X();
                vertex.Y(PointLayer.GetFeature(i).GetFieldAsDouble("POINT_Y"));
                if (ymin > vertex.Y())
                    ymin = vertex.Y();
                if (ymax < vertex.Y())
                    ymax = vertex.Y();
                vertex.Z(DEMServer.GetElevation(vertex.X(), vertex.Y()));
                if (0 != PointLayer.GetFeature(i).GetFieldAsInteger("inclinatio"))
                    vertex.occurrence.inclination = (double) PointLayer.GetFeature(i).GetFieldAsInteger("inclinatio");
                vertex.type = 0;

                //将产状点产状赋值给相同地层最近边界点
                //这里先随机给几个边界点添加倾角信息，模拟这个赋值过程，之后插值出其他边界点倾角
                if (0 == i % ((int)PointLayer.GetFeatureCount(0) / 3) && 0 == si)
                    vertex.occurrence.dip = Math.PI / 6;
                //vertex.occurrence.dip = (25 + rd.Next() % (40 - 25)) / 180.0 * Math.PI;
                else if (0 == i % ((int)PointLayer.GetFeatureCount(0) / 3) && 1 == si)
                    vertex.occurrence.dip = Math.PI / 3;
                    //vertex.occurrence.dip = (45 + rd.Next() % (55 - 45)) / 180.0 * Math.PI;

                vcUp.addVer(vertex);   
            }
            if (vcUp.Count < 3)
                MessageBox.Show("地层产状点不足！");
            fornix.setRange(xmin, ymin, xmax, ymax);
        }

        /// <summary>
        /// 计算地层边界点角度与曲率
        /// </summary>
        /// <param name="fornix"></param>
        private static void calVertexAC(Fornix fornix)
        {
            Vertex prev, curv, nextv;
            for (int i = 1; i < fornix.outSide.upvers.Count - 1; ++i)
            {
                prev = fornix.outSide.getUpver(i - 1);
                curv = fornix.outSide.getUpver(i);
                nextv = fornix.outSide.getUpver(i + 1);
                curv.angle = curv.calAngle(prev, nextv) > 0.0 ? Math.Acos(curv.calAngle(prev, nextv)) : Math.PI - Math.Acos(curv.calAngle(prev, nextv));
                curv.calCurvature(prev, nextv);
            }

            prev = fornix.outSide.getUpver(fornix.outSide.upvers.Count - 1);
            curv = fornix.outSide.getUpver(0);
            nextv = fornix.outSide.getUpver(1);
            curv.angle = curv.calAngle(prev, nextv) > 0.0 ? Math.Acos(curv.calAngle(prev, nextv)) : Math.PI - Math.Acos(curv.calAngle(prev, nextv));
            curv.calCurvature(prev, nextv);

            prev = fornix.outSide.getUpver(fornix.outSide.upvers.Count - 2);
            curv = fornix.outSide.getUpver(fornix.outSide.upvers.Count - 1);
            nextv = fornix.outSide.getUpver(0);
            curv.angle = curv.calAngle(prev, nextv) > 0.0 ? Math.Acos(curv.calAngle(prev, nextv)) : Math.PI - Math.Acos(curv.calAngle(prev, nextv));
            curv.calCurvature(prev, nextv);

            ObjWriter.writeAC(fornix);
        }

        /// <summary>
        /// 计算地层边界点产状
        /// </summary>
        /// <param name="fornix"></param>
        private static void calOccurrence(Fornix fornix)
        {
            //求地层线等高线交点，并据此计算一部分边界点倾向和倾角
            StrataContourIntersect(fornix);
            //根据已知产状点和地层等高线交点的倾角值，线性插值计算其他边界点倾角
            calDips(fornix);
            //先行计算地层界线噪声部位的倾向
            calNoiseInc(fornix);
            //根据地层等高线交点倾向值和噪声部位倾向值，线性插值剩下的边界点倾向
            calInclinations(fornix);
            //后处理边界点倾向，若上述计算导致其倾向朝内，则改用外接圆法计算倾向
            calInclinationsAfter(fornix);
        }

        /// <summary>
        /// 求地层界线与等高线的交点，并用三点法求交点的倾向值，再把倾向值赋予其前一个地层边界点
        /// </summary>
        /// <param name="fornix"></param>
        private static void StrataContourIntersect(Fornix fornix)
        {
            VertexCollection vc;
            Contour contour;
            InsectVer insection;
            List<InsectVer> insections = new List<InsectVer>();
            int i, j, k;
            for (i = 0; i < _contours.Count;++i )
            {
                contour = _contours[i];
                if(fornix.fastJudge(contour))
                    continue;
                vc = fornix.outSide.upvers;
                for (j = 0; j < contour.Count - 1; ++j)
                {
                    for (k = 0; k < vc.Count - 1; ++k)
                    {
                        if (Tools.InsectionJudge(contour.getVer(j), contour.getVer(j + 1), vc.getVer(k), vc.getVer(k + 1)))
                        {
                            insection = Tools.GetCrossPoint(contour.getVer(j), contour.getVer(j + 1), vc.getVer(k), vc.getVer(k + 1));
                            if (null == insection)
                                continue;
                            //交点记录邻接边界点ID
                            insection.VerIdx = k;
                            //交点同时放入交点集合中
                            insections.Add(insection);
                        }
                    }
                }
            }
            //求交点集合中各个交点的倾向值
            calInsectInc(insections, fornix);
        }

        /// <summary>
        /// 三点法求取地层等高线交点倾向
        /// </summary>
        /// <param name="insections"></param>
        private static void calInsectInc(List<InsectVer> insections, Fornix fornix)
        {
            if (3 > insections.Count)
            {
                MessageBox.Show("当前地层与等高线交点数量不足三个！");
                return;
            }
            InsectVer v0, v1, v2;
            int i;
            for (i = 0; i < insections.Count - 2; ++i)
            {
                v0 = insections[i];
                v1 = insections[i + 1]; 
                v2 = insections[i + 2];
                if (!(v0.Z() == v1.Z() && v0.Z() == v1.Z()))
                {
                    //计算交点产状
                    v1.calOccuurence1(v0, v2);
                    //将交点产状赋值给其上邻接边界点
                    fornix.outSide.upvers.getVer(v1.VerIdx).occurrence.inclination = v1.occurrence.inclination;
                    if (v1.occurrence.dip > Math.PI / 6 && v1.occurrence.dip < Math.PI / 12 * 5)//避免个别倾角过小，导致侧面过长畸形
                        fornix.outSide.upvers.getVer(v1.VerIdx).occurrence.dip = v1.occurrence.dip;
                }
            }
            v0 = insections[insections.Count - 2];
            v1 = insections[insections.Count - 1];
            v2 = insections[0];
            if (!(v0.Z() == v1.Z() && v0.Z() == v1.Z()))
            {
                v1.calOccuurence1(v0, v2);
                fornix.outSide.upvers.getVer(v1.VerIdx).occurrence.inclination = v1.occurrence.inclination;
                if (v1.occurrence.dip > Math.PI / 6)
                    fornix.outSide.upvers.getVer(v1.VerIdx).occurrence.dip = v1.occurrence.dip;
            }
            v0 = insections[insections.Count - 1];
            v1 = insections[0];
            v2 = insections[1];
            if (!(v0.Z() == v1.Z() && v0.Z() == v1.Z()))
            {
                v1.calOccuurence1(v0, v2);
                fornix.outSide.upvers.getVer(v1.VerIdx).occurrence.inclination = v1.occurrence.inclination;
                if (v1.occurrence.dip > Math.PI / 6)
                    fornix.outSide.upvers.getVer(v1.VerIdx).occurrence.dip = v1.occurrence.dip;
            }

        }

        /// <summary>
        /// 生成穹窿模型侧面
        /// </summary>
        /// <param name="fornix"></param>
        /// <param name="preFornix"></param>
        private static void createSide(Fornix fornix, Fornix preFornix)
        {
            //根据产状和贝塞尔系数生成几圈侧边点
            _verNum += fornix.createOutSideCircles(_downElevation, Bessel_rad_cof, Bessel_num_control, Bessel_num_collected);
            //根据产状生成底面顶点
            //fornix.createOutSideLowvers(createDownVer(fornix));//底面顶点集合
            fornix.createOutSideLowversFromSide();
            //一层层生成地层外侧面
            _verNum -= fornix.createOutSideByBessel(Bessel_num_collected);
            //_verNum -= fornix.createOutSide();
            //生成前一地层的内侧面
            if (preFornix != null)
                preFornix.createInSide(fornix);
        }

        /// <summary>
        /// 由倾角已知的边界点插值计算其他边界点倾角值
        /// </summary>
        /// <param name="fornix"></param>
        private static void calDips(Fornix fornix)
        {
            Vertex vertexUp;
            int preIdx = -1, firstIdx = -1;
            double preDip = -1.0, curDip = -1.0, firstDip = -1.0;
            for (int i = 0; i < fornix.outSide.countUpVers(); ++i)
            {
                vertexUp = fornix.outSide.getUpver(i);
                if (0.0 != vertexUp.occurrence.dip)
                {
                    curDip = vertexUp.occurrence.dip;
                    if (double.IsNaN(curDip))
                        Console.WriteLine("");
                    if (-1 != preIdx)
                        for (int j = preIdx + 1; j < i; ++j)
                            fornix.outSide.getUpver(j).occurrence.dip = preDip + (curDip - preDip) / (i - preIdx) * (j - preIdx);//线性插值
                    else
                    {
                        firstIdx = i;
                        firstDip = curDip;
                    }
                    preIdx = i;
                    preDip = curDip;
                }
            }
            //处理最后一个已知点和第一个已知点中间的产状点倾角
            for (int i = preIdx + 1; i < fornix.outSide.countUpVers(); ++i)
                fornix.outSide.getUpver(i).occurrence.dip = preDip + (firstDip - preDip) / (fornix.outSide.countUpVers() - preIdx + firstIdx - 1) * (i - preIdx);
            for (int i = 0; i < firstIdx; ++i)
                fornix.outSide.getUpver(i).occurrence.dip = preDip + (firstDip - preDip) / (fornix.outSide.countUpVers() - preIdx + firstIdx - 1) * (fornix.outSide.countUpVers() - preIdx + i);

        }

        /// <summary>
        /// 计算地层噪声部位的倾向
        /// </summary>
        /// <param name="fornix"></param>
        private static void calNoiseInc(Fornix fornix)
        {
            
        }

        /// <summary>
        /// 由倾向已知的边界点插值计算其他边界点倾向值
        /// </summary>
        /// <param name="fornix"></param>
        private static void calInclinations(Fornix fornix)
        {
            Vertex vertexUp;
            int preIdx = -1, firstIdx = -1;
            double preInc = -1.0, curInc = -1.0, firstInc = -1.0;
            for (int i = 0; i < fornix.outSide.countUpVers(); ++i)
            {
                vertexUp = fornix.outSide.getUpver(i);
                if (0.0 != vertexUp.occurrence.inclination)
                {
                    curInc = vertexUp.occurrence.inclination;
                    if (-1 != preIdx)
                        for (int j = preIdx + 1; j < i; ++j)
                            fornix.outSide.getUpver(j).occurrence.inclination = preInc + (curInc - preInc) / (i - preIdx) * (j - preIdx);//线性插值
                    else
                    {
                        firstIdx = i;
                        firstInc = curInc;
                    }
                    preIdx = i;
                    preInc = curInc;
                }
            }
            //处理最后一个已知点和第一个已知点中间的产状点倾向
            for (int i = preIdx + 1; i < fornix.outSide.countUpVers(); ++i)
                fornix.outSide.getUpver(i).occurrence.inclination = preInc + (firstInc - preInc) / (fornix.outSide.countUpVers() - preIdx + firstIdx - 1) * (i - preIdx);
            for (int i = 0; i < firstIdx; ++i)
                fornix.outSide.getUpver(i).occurrence.inclination = preInc + (firstInc - preInc) / (fornix.outSide.countUpVers() - preIdx + firstIdx - 1) * (fornix.outSide.countUpVers() - preIdx + i);

        }

        /// <summary>
        /// 后处理边界点倾向，若上述计算导致其倾向朝内，则改用外接圆法计算倾向
        /// </summary>
        /// <param name="fornix"></param>
        private static void calInclinationsAfter(Fornix fornix)
        {
            Vertex prev, curv, nextv, incv;
            for (int i = 1; i < fornix.outSide.upvers.Count - 1; ++i)
            {
                prev = fornix.outSide.getUpver(i - 1);
                curv = fornix.outSide.getUpver(i);
                nextv = fornix.outSide.getUpver(i + 1);
                incv = new Vertex();
                incv.X(curv.X() + Math.Sin(curv.occurrence.inclination));
                incv.Y(curv.Y() + Math.Cos(curv.occurrence.inclination));
                if (0.0 > curv.calVector(prev, nextv) && (0.0 > curv.calVector(incv, nextv) || 0.0 > curv.calVector(prev, incv)))
                    curv.calOccuurence(prev, nextv);
                else if (0.0 < curv.calVector(prev, nextv) && (0.0 > curv.calVector(incv, nextv) || 0.0 < curv.calVector(incv, prev)))
                    curv.calOccuurence(prev, nextv);
            }

            prev = fornix.outSide.getUpver(fornix.outSide.upvers.Count - 1);
            curv = fornix.outSide.getUpver(0);
            nextv = fornix.outSide.getUpver(1);
            incv = new Vertex();
            incv.X(curv.X() + Math.Sin(curv.occurrence.inclination));
            incv.Y(curv.Y() + Math.Cos(curv.occurrence.inclination));
            if (0.0 > curv.calVector(prev, nextv) && (0.0 > curv.calVector(incv, nextv) || 0.0 > curv.calVector(prev, incv)))
                curv.calOccuurence(prev, nextv);
            else if (0.0 < curv.calVector(prev, nextv) && (0.0 > curv.calVector(incv, nextv) || 0.0 < curv.calVector(incv, prev)))
                curv.calOccuurence(prev, nextv);

            prev = fornix.outSide.getUpver(fornix.outSide.upvers.Count - 2);
            curv = fornix.outSide.getUpver(fornix.outSide.upvers.Count - 1);
            nextv = fornix.outSide.getUpver(0);
            incv = new Vertex();
            incv.X(curv.X() + Math.Sin(curv.occurrence.inclination));
            incv.Y(curv.Y() + Math.Cos(curv.occurrence.inclination));
            if (0.0 > curv.calVector(prev, nextv) && (0.0 > curv.calVector(incv, nextv) || 0.0 > curv.calVector(prev, incv)))
                curv.calOccuurence(prev, nextv);
            else if (0.0 < curv.calVector(prev, nextv) && (0.0 > curv.calVector(incv, nextv) || 0.0 < curv.calVector(incv, prev)))
                curv.calOccuurence(prev, nextv);
        }

        /// <summary>
        /// 由边界点产状生成侧边线控制点
        /// </summary>
        /// <param name="fornix"></param>
        /// <param name="vcDown"></param>
        private static VertexCollection createLowerVer(Fornix fornix)
        {
            VertexCollection vcDown = new VertexCollection();
            Vertex vertexUp, vertexDown;
            for (int i = 0; i < fornix.outSide.countUpVers(); ++i)
            {
                vertexUp = fornix.outSide.getUpver(i);
                vertexDown = new Vertex();
                vertexDown.createDownVer(vertexUp, _downElevation);
                vertexDown.ID = _verNum++;
                vertexDown.corIdx = i;
                vertexUp.corIdx = i;
                vcDown.addVer(vertexDown);
            }
            return vcDown;
        }

        /// <summary>
        /// 由边界点产状计算直接生成底边顶点
        /// </summary>
        /// <param name="fornix"></param>
        /// <param name="vcDown"></param>
        private static VertexCollection createDownVer(Fornix fornix)
        {
            VertexCollection vcDown = new VertexCollection();
            Vertex vertexUp, vertexDown;
            for (int i = 0; i < fornix.outSide.countUpVers(); ++i)
            {
                vertexUp = fornix.outSide.getUpver(i);
                vertexDown = new Vertex();
                vertexDown.createDownVer(vertexUp, _downElevation);
                vertexDown.ID = _verNum++;
                vertexDown.corIdx = i;
                vertexUp.corIdx = i;
                vcDown.addVer(vertexDown);
            }
            return vcDown;
        }

        /// <summary>
        /// 读入SHP文件
        /// </summary>
        /// <param name="shpFilePath"></param>
        private static void ReadShp(int status, string shpFilePath)
        {
            _ds = Ogr.Open(shpFilePath, 1);//0表示只读，1表示可修改  
            if (_ds == null) { MessageBox.Show("打开文件【{0}】失败！", shpFilePath); return; }
            // 获取第一个图层
            Layer curLayer = _ds.GetLayerByIndex(0);
            if (curLayer == null) { MessageBox.Show("获取第{0}个图层失败！ n", "0"); return; }
            if (0 == status)
                _layers.Add(curLayer);
            else if (1 == status)
                _contourLayer = curLayer;
            
        }

        /// <summary>
        /// 输出模型底面shp
        /// </summary>
        /// <param name="fornixs"></param>
        /// <param name="shppath"></param>
        public static void writeDownSHP(List<Fornix> fornixs, string shppath)
        {
            string strVectorFile;
            // 为了支持中文路径，请添加下面这句代码
            OSGeo.GDAL.Gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
            // 为了使属性表字段支持中文，请添加下面这句
            OSGeo.GDAL.Gdal.SetConfigOption("SHAPE_ENCODING", "");
            // 注册所有的驱动
            Ogr.RegisterAll();
            for (int i = 0; i < fornixs.Count; ++i)
            {
                Fornix fornix = fornixs[i];
                strVectorFile = shppath + fornix.name + "_down.shp";
                //创建数据，创建ESRI的shp文件
                string strDriverName = "ESRI Shapefile";
                Driver oDriver = Ogr.GetDriverByName(strDriverName);
                if (oDriver == null)
                {
                    Console.WriteLine("%s 驱动不可用！\n", strVectorFile);
                    return;
                }

                // 步骤1、创建数据源
                DataSource oDS = oDriver.CreateDataSource(strVectorFile, null);
                if (oDS == null)
                {
                    Console.WriteLine("创建矢量文件【%s】失败！", strVectorFile);
                    return;
                }
                //步骤2、创建空间坐标系
                OSGeo.OSR.SpatialReference oSRS = new OSGeo.OSR.SpatialReference("");
                //oSRS.SetWellKnownGeogCS("WGS84");
                //步骤3、创建图层，并添加坐标系，创建一个多边形图层(wkbGeometryType.wkbUnknown,存放任意几何特征)
                Layer oLayer = oDS.CreateLayer("TestPolygon", oSRS, wkbGeometryType.wkbUnknown, null);
                if (oLayer == null)
                {
                    Console.WriteLine("图层创建失败！");
                    return;
                }

                // 步骤4、下面创建属性表
                FieldDefn oFieldPlotArea = new FieldDefn("PlotArea", FieldType.OFTString);          // 先创建一个叫PlotArea的属性
                oFieldPlotArea.SetWidth(100);
                // 步骤5、将创建的属性表添加到图层中
                oLayer.CreateField(oFieldPlotArea, 1);
                //步骤6、定义一个特征要素oFeature(特征要素包含两个方面1.属性特征2.几何特征)
                FeatureDefn oDefn = oLayer.GetLayerDefn();
                Feature oFeature = new Feature(oDefn);    //建立了一个特征要素并将指向图层oLayer的属性表
                //步骤7、设置属性特征的值
                //oFeature.SetField(0, area.ToString());
                oFeature.SetField(0, "1234");

                string wkt = "POLYGON((";
                for (int j = 0; j < fornix.outSide.lowvers.Count; ++j)
                {
                    wkt += fornix.outSide.getDownver(j).X().ToString() + " " + fornix.outSide.getDownver(j).Y().ToString();
                    if (j != fornix.outSide.lowvers.Count - 1)
                        wkt += ",";
                }
                wkt += "))";
                OSGeo.OGR.Geometry geomTriangle = OSGeo.OGR.Geometry.CreateFromWkt(wkt);//创建一个几何特征
                //步骤8、设置几何特征
                oFeature.SetGeometry(geomTriangle);
                //步骤9、将特征要素添加到图层中
                oLayer.CreateFeature(oFeature);
                //Debug.WriteLine("数据集创建完成！");
            }

            
        }


        private static void SetShp(Layer PointLayer, Fornix fornix)
        {
            if (-1 == PointLayer.FindFieldIndex("trend", 0))
            {
                FieldDefn oFieldName0 = new FieldDefn("trend", FieldType.OFTReal);
                oFieldName0.SetWidth(50);
                oFieldName0.SetPrecision(7);
                PointLayer.CreateField(oFieldName0, 1);
            }
            if (-1 == PointLayer.FindFieldIndex("incli", 0))
            {
                FieldDefn oFieldName1 = new FieldDefn("incli", FieldType.OFTReal);
                oFieldName1.SetWidth(50);
                oFieldName1.SetPrecision(7);
                PointLayer.CreateField(oFieldName1, 1);
            }
            if (-1 == PointLayer.FindFieldIndex("dip", 0))
            {
                FieldDefn oFieldName2 = new FieldDefn("dip", FieldType.OFTReal);
                oFieldName2.SetWidth(50);
                oFieldName2.SetPrecision(7);
                PointLayer.CreateField(oFieldName2, 1);
            }

            for (int i = 0; i < (int)PointLayer.GetFeatureCount(0); ++i)
            {
                Vertex vertex = fornix.outSide.upvers.getVer(i);
                Feature pointFeature = PointLayer.GetFeature(i);
                pointFeature.SetField("trend", vertex.occurrence.trend);
                pointFeature.SetField("incli", vertex.occurrence.inclination);
                pointFeature.SetField("dip", vertex.occurrence.dip);
                PointLayer.SetFeature(pointFeature);//更改其值
                pointFeature.Dispose();//释放对象
            }
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
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        /// <returns>交点</returns>
        private static InsectVer GetCrossPoint(Vertex a, Vertex b, Vertex c, Vertex d)
        {
            InsectVer crossPoint = new InsectVer();
            double A = (b.X() - a.X()) / (b.Y() - a.Y());
            double B = (d.X() - c.X()) / (d.Y() - c.Y());
            double y = (c.X() - a.X() + A * a.Y() - B * c.Y()) / (A - B);
            double x = A * y - A * a.Y() + a.X();
            if (double.IsNaN(x) || double.IsNaN(y))
                return null;
            crossPoint.X(x);
            crossPoint.Y(y);
            crossPoint.Z((a.Z() + b.Z() + c.Z() + d.Z()) / 4);
            crossPoint.ID = b.ID;
            return crossPoint;
        }

    }
}
