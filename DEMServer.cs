using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Threading.Tasks;
using OSGeo.GDAL;

namespace FornixModelingGDAL
{
    class DEMServer
    {
        private static Dataset _ds;
        private static Band _band;
        private static double[] _adfGeoTransform = new double[6];
        private static int srcWidth = -1;
        private static int srcHeight = -1;
        private static int bandCount = -1;
        private static double[] lim = new double[2];

        public static VertexCollection getVersFromDEM(string imgpath)
        {

            VertexCollection vc = new VertexCollection();
            //栅格范围：left，top，right，bottom
            _ds = Gdal.Open(imgpath, Access.GA_ReadOnly);
            srcWidth = _ds.RasterXSize;
            srcHeight = _ds.RasterYSize;
            bandCount = _ds.RasterCount;
            int[] bandArray = new int[bandCount];
            for (int i = 0; i < bandCount; i++)
                bandArray[i] = i + 1;
            double[] dataArray = new double[srcWidth * srcHeight * bandCount];
            double[] dataArray1 = new double[srcWidth * srcHeight * bandCount];
            double x, y;
            _ds.ReadRaster(0, 0, srcWidth, srcHeight, dataArray, srcWidth, srcHeight, bandCount, bandArray, 0, 0, 0);
            //RasterWriter.writeFornixObj(dataArray, srcWidth, srcHeight, @"E:\Users\LiuXianyu\Documents\ExperimentData\myProject\FornixModelingGDAL\Data\LingYanShan\Export\Raster.txt");
            //获取坐标变换系数
            //0左上角x坐标
            //1东西方向分辨率
            //2旋转角度, 0表示图像 "北方朝上"
            //3左上角y坐标
            //4旋转角度, 0表示图像 "北方朝上"
            //5南北方向分辨率
            _ds.GetGeoTransform(_adfGeoTransform);
            Band band = _ds.GetRasterBand(1);
            band.ComputeRasterMinMax(lim, 0);

            //VIP法筛选DEM点
            //getInnerPointsByVIP(dataArray, dataArray1);

            //固定步长法筛选DEM点
            getInnerPointsByStep(dataArray, dataArray1);

            //根据矩阵数据生成点集
            for (int i = 0; i < srcHeight; ++i)
                for (int j = 0; j < srcWidth; ++j)
                {
                    if (-1.0 != dataArray1[i * srcWidth + j])
                    {
                        Vertex ver = new Vertex();
                        getCoordinateFromMatrix(dataArray1, i, j, out x, out y);
                        ver.X(x);
                        ver.Y(y);
                        //ver.Z(dataArray1[i * srcWidth + j]);
                        ver.Z(GetElevation(x, y));
                        ver.innerPoint = true;
                        vc.addVer(ver);
                    }
                }
            return vc;
        }

        private static double calAgvH(double[] dataArray, int i, int j, int srcWidth)
        {
            double[] elevation = new double[9];
            double[] height = new double[4];
            double resolution = _adfGeoTransform[1];
            int idx = 0;
            for (int s = i - 1; s <= i + 1; ++s)
                for (int t = j - 1; t <= j + 1; ++t)
                    elevation[idx++] = dataArray[s * srcWidth + t];
            height[0] = (resolution * 2 * elevation[4] - resolution * elevation[0] - resolution * elevation[8]) * Math.Sqrt(2)
                / Math.Sqrt(Math.Pow(resolution * 2 * Math.Sqrt(2), 2) + Math.Pow(elevation[0] - elevation[8], 2));
            height[1] = (resolution * 2 * elevation[4] - resolution * elevation[1] - resolution * elevation[7])
                / Math.Sqrt(Math.Pow(resolution * 2, 2) + Math.Pow(elevation[1] - elevation[7], 2));
            height[2] = (resolution * 2 * elevation[4] - resolution * elevation[2] - resolution * elevation[6]) * Math.Sqrt(2)
                / Math.Sqrt(Math.Pow(resolution * 2 * Math.Sqrt(2), 2) + Math.Pow(elevation[2] - elevation[6], 2));
            height[3] = (resolution * 2 * elevation[4] - resolution * elevation[3] - resolution * elevation[5])
                / Math.Sqrt(Math.Pow(resolution * 2, 2) + Math.Pow(elevation[3] - elevation[5], 2));
            return (height[0] + height[1] + height[2] + height[3]) / 4;
        }

        private static void getCoordinateFromMatrix(double[] dataArray, int i, int j, out double x, out double y)
        {
            x = _adfGeoTransform[0] + _adfGeoTransform[1] * (i + 0.5);
            y = _adfGeoTransform[3] + _adfGeoTransform[5] * (j + 0.5);
        }

        /// <summary>
        /// VIP法筛选DEM点
        /// </summary>
        /// <param name="dataArray"></param>
        /// <param name="dataArray1"></param>
        private static void getInnerPointsByVIP(double[] dataArray, double[] dataArray1)
        {
            double avgh, total_avgh = 0.0, avg_avgh;
            int num = 0;
            for (int i = 0; i < srcWidth; ++i)
                for (int j = 0; j < srcHeight; ++j)
                {
                    if (i == 0 || i == srcWidth - 1 || j == 0 || j == srcHeight - 1)
                    {
                        dataArray1[i * srcWidth + j] = -1.0;
                        continue;
                    }
                    avgh = calAgvH(dataArray, i, j, srcWidth);
                    total_avgh += Math.Abs(avgh);
                    if (dataArray[i * srcWidth + j] > 23 && (avgh > 1.5 || avgh < -1.5))//大于阈值，保留
                    {
                        dataArray1[i * srcWidth + j] = dataArray[i * srcWidth + j];
                        ++num;
                    }
                    else
                        dataArray1[i * srcWidth + j] = -1.0;
                }
            avg_avgh = total_avgh / (srcWidth - 1) / (srcHeight - 1);
        }

        /// <summary>
        /// 固定步长法筛选DEM点
        /// </summary>
        /// <param name="dataArray"></param>
        /// <param name="dataArray1"></param>
        private static void getInnerPointsByStep(double[] dataArray, double[] dataArray1)
        {
            double total_avgh = 0.0, avg_avgh;
            int stepWidth = (int)Math.Round((double)srcWidth / 50);
            int stepHeight = (int)Math.Round((double)srcHeight / 50);

            for (int i = 0; i < srcHeight; ++i)
                for (int j = 0; j < srcWidth; ++j)
                {
                    if (!(i % stepWidth == 0 && j % stepHeight == 0))
                    {
                        if (dataArray[i * srcWidth + j] > lim[1] - 0.01)
                            dataArray1[i * srcWidth + j] = dataArray[i * srcWidth + j];
                        else
                            dataArray1[i * srcWidth + j] = -1.0;
                        continue;
                    }
                    dataArray1[i * srcWidth + j] = dataArray[i * srcWidth + j];
                    total_avgh += dataArray[i * srcWidth + j];
                }
            avg_avgh = total_avgh / (srcWidth - 1) / (srcHeight - 1);
        }


        /// <summary>
        /// 获取DEM上指定点的高程值
        /// </summary>
        /// <param name="dProjX"></param>
        /// <param name="dProjY"></param>
        /// <returns></returns>
        public static double GetElevation(double dProjX, double dProjY)
        {
            try
            {
                _band = _ds.GetRasterBand(1);
                //获取图像的尺寸               
                int width = _band.XSize;
                int height = _band.YSize;

                //获取行列号
                double dTemp = _adfGeoTransform[1] * _adfGeoTransform[5] - _adfGeoTransform[2] * _adfGeoTransform[4];
                double dCol = 0.0, dRow = 0.0;
                dCol = (_adfGeoTransform[5] * (dProjX - _adfGeoTransform[0]) -
                    _adfGeoTransform[2] * (dProjY - _adfGeoTransform[3])) / dTemp + 0.5;
                dRow = (_adfGeoTransform[1] * (dProjY - _adfGeoTransform[3]) -
                    _adfGeoTransform[4] * (dProjX - _adfGeoTransform[0])) / dTemp + 0.5;
                int dc = Convert.ToInt32(dCol);
                int dr = Convert.ToInt32(dRow);


                //获取DEM数值到一维数组
                double[] data = new double[1 * 1];
                CPLErr err = _band.ReadRaster(dc, dr, 1, 1, data, 1, 1, 0, 0);
                _band.Dispose();
                double elvate = data[0];
                return elvate;
            }
            catch
            {
                return 0.0;
            }
        }

        /// <summary>
        /// 提取等高线
        /// </summary>
        /// <param name="shpPath">输出路径</param>
        /// <param name="dfContourInterval">等高线间隔</param>
        /// <returns></returns>
        //bool CalContourByGdal(string shpPath, double dfContourInterval)
        //{
            //if(_ds == null)  
            //{  
            //    MessageBox.Show("读取图像失败！");  
            //    return false;  
            //}

            ////读取图像数据波段
            //_band = _ds.GetRasterBand(1);
            //float *pData = new float[srcWidth*srcHeight]();

            ////获取行列号
            //double dTemp = _adfGeoTransform[1] * _adfGeoTransform[5] - _adfGeoTransform[2] * _adfGeoTransform[4];
            //double dCol = 0.0, dRow = 0.0;
            //dCol = (_adfGeoTransform[5] * (dProjX - _adfGeoTransform[0]) -
            //    _adfGeoTransform[2] * (dProjY - _adfGeoTransform[3])) / dTemp + 0.5;
            //dRow = (_adfGeoTransform[1] * (dProjY - _adfGeoTransform[3]) -
            //    _adfGeoTransform[4] * (dProjX - _adfGeoTransform[0])) / dTemp + 0.5;
            //int dc = Convert.ToInt32(dCol);
            //int dr = Convert.ToInt32(dRow);


            ////获取DEM数值到一维数组
            //double[] data = new double[1 * 1];
            //CPLErr err = _band.ReadRaster(dc, dr, 1, 1, data, 1, 1, 0, 0);

            //CPLErr err = _band->RasterIO(GF_Read, 0, 0, srcWidth, srcHeight, pData, srcWidth, srcHeight, GDT_Float32, 0, 0);
            //if(err == CPLErr.CE_Failure)
            //{
            //    MessageBox.Show("读取DEM图像数据时出错！");
            //    return false;
            //}

            ////判断图像中是否有异常值，并获取异常值实际值
            //float fNoData = 0;
            //int nIdx;
            //for(int i=0; i<srcHeight; i++) 
            //{
            //    for(int j=0; j<srcWidth; j++)
            //    {
            //        nIdx = i*srcWidth + j;
            //        if(pData[nIdx] <= -9999)
            //        {
            //            fNoData = pData[nIdx];
            //        }
            //    }
            //}


            ////创建矢量图
            //const char *pszDriverName = "ESRI Shapefile";
            //GDALDriver *poDriver;
            //GDALAllRegister();
            //poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
            //if( poDriver == NULL )
            //{
            //    printf( "%s driver not available.\n", pszDriverName );
            //    exit( 1 );
            //}
            //GDALDataset *poDS;
            //poDS = poDriver->Create( shpPath.c_str(), 0, 0, 0, GDT_Unknown, NULL ); //创建数据源
            //if( poDS == NULL )
            //{
            //    printf( "Creation of output file failed.\n" );
            //    exit( 1 );
            //}
            //OGRLayer *poLayer;
            //OGRSpatialReference *poSpatialRef = new OGRSpatialReference(_ds->GetProjectionRef());
            //poLayer = poDS->CreateLayer( "Elevation", poSpatialRef, wkbLineString, NULL ); //创建图层
            //if( poLayer == NULL )
            //{
            //    printf( "Layer creation failed.\n" );
            //    exit( 1 );
            //}


            //OGRFieldDefn ofieldDef1("Elevation", OFTInteger);    //在矢量图中创建高程值字段
            //if(poLayer->CreateField(&ofieldDef1) != OGRERR_NONE)
            //{
            //    cout<<"创建矢量图层属性表失败！"<<endl;
            //    GDALClose((GDALDatasetH)poDS);
            //    GDALClose(_ds);
            //    return 0;
            //}

            ////根据图像波段生成矢量图等高线
            //if(fNoData == 0) 
            //    GDALContourGenerate(_band, dfContourInterval, 0, 0, NULL, false, 0, (OGRLayerH)poLayer, -1, 0, NULL, NULL);
            //else //有异常值时，不对异常值进行处理
            //    GDALContourGenerate(_band, dfContourInterval, 0, 0, NULL, true, fNoData, (OGRLayerH)poLayer, -1, 0, NULL, NULL);

            //GDALClose(poDS);
            //GDALClose(_ds);


            //if(pData != NULL)
            //{
            //    delete[] pData;
            //    pData = NULL;
            //}

        //}

        /// <summary>
        /// GDAL方法提取等值线
        /// </summary>
        /// <param name="imgPath"></param>
        /// <param name="shpPath"></param>
        /// <param name="dfContourInterval"></param>等高线间隔
        /// <returns></returns>
        //public void CalContourByGdal(string imgPath, string shpPath, string fileName, double dfContourInterval)
        //{
        //    Dataset pInDataset = Gdal.Open(imgPath, Access.GA_ReadOnly);

        //    if (pInDataset == null)
        //    {
        //        Console.WriteLine("读取图像失败！");
        //        //return false;
        //    }

        //    int nDemWidth = pInDataset.RasterXSize; //获取图像宽  
        //    int nDemHeight = pInDataset.RasterYSize; //获取图像高  
        //    int Count = pInDataset.RasterCount; //波段数  


        //    //读取图像数据波段
        //    Band pInRasterBand = pInDataset.GetRasterBand(1);
        //    float[] pData = new float[nDemWidth * nDemHeight];


        //    CPLErr err = pInRasterBand.ReadRaster(0, 0, nDemWidth, nDemHeight, pData, nDemWidth, nDemHeight, 0, 0);
        //    if (err == CPLErr.CE_Failure)
        //    {
        //        Console.WriteLine("读取DEM图像数据时出错！");
        //        Gdal.GDALDestroyDriverManager();
        //        //return false;
        //    }
        //    //判断图像中是否有异常值，并获取异常值实际值
        //    float fNoData = 0;
        //    int nIdx;
        //    for (int i = 0; i < nDemHeight; i++)
        //    {
        //        for (int j = 0; j < nDemWidth; j++)
        //        {
        //            nIdx = i * nDemWidth + j;
        //            if (pData[nIdx] <= -9999)
        //            {
        //                fNoData = pData[nIdx];
        //            }
        //        }
        //    }

        //    //创建矢量图           
        //    //OSGeo.GDAL.Driver poDriver;
        //    //Dataset poDS;           
        //    //poDS = poDriver.Create(shpPath, 0, 0, 0, DataType.GDT_Unknown, null); //创建数据源
        //    OSGeo.OGR.DataSource poDS;
        //    poDS = oDerive.CreateDataSource(shpPath + "\\" + fileName + ".shp", null);
        //    if (poDS == null)
        //    {
        //        Console.WriteLine("Creation of output file failed.\n");
        //        //exit(1);System.Environment.Exit(-1);
        //        //return false;
        //    }
        //    SpatialReference poSpatialRef = new SpatialReference(pInDataset.GetProjectionRef());
        //    oLayer = poDS.CreateLayer(fileName, poSpatialRef, wkbGeometryType.wkbLineString, null); //创建图层            
        //    if (oLayer == null)
        //    {
        //        Console.WriteLine("Layer creation failed.\n");
        //        //exit(1);
        //        //return false;
        //    }

        //    FieldDefn ofieldDef = new FieldDefn("Elevation", FieldType.OFTInteger);    //在矢量图中创建高程值字段
        //    if (oLayer.CreateField(ofieldDef, 1) != Ogr.OGRERR_NONE)
        //    {
        //        Console.WriteLine("创建矢量图层属性表失败！");
        //        //GDALClose((GDALDatasetH)poDS);
        //        poDS.Dispose();
        //        //GDALClose(pInDataset);
        //        pInDataset.Dispose();
        //        //return true;
        //    }

        //    //根据图像波段生成矢量图等高线
        //    if (fNoData == 0)
        //        Gdal.ContourGenerate(pInRasterBand, dfContourInterval, 0, 0, null, -1, 0, oLayer, -1, 0, null, null);//*int bUseNoData 是否使用无效值 false
        //    else //有异常值时，不对异常值进行处理
        //        Gdal.ContourGenerate(pInRasterBand, dfContourInterval, 0, 0, null, 1, fNoData, oLayer, -1, 0, null, null);//true

        //    //GDALClose(poDS);
        //    poDS.Dispose();
        //    pInDataset.Dispose();
        //    //GDALClose(pInDataset);

        //    if (pData != null)
        //    {
        //        Array.Clear(pData, 0, pData.Length);
        //        pData = null;
        //    }

        //}

    }
}
