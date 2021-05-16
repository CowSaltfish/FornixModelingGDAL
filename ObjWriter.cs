using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace FornixModelingGDAL
{
    class ObjWriter
    {
        /// <summary>
        /// 导出穹窿模型
        /// </summary>
        /// <param name="fornixs"></param>
        /// <param name="path"></param>
        /// <param name="verNum">边界点个数</param>
        public static void writeFornixObj(List<Fornix> fornixs, string path,int verNum)
        {
            StreamWriter sw = File.CreateText(path);
            int i;
            foreach (Fornix fornix in fornixs)
            {
                //生成顶面外顶点
                for (i = 0; i < fornix.outSide.upvers.Count; ++i)
                {
                    sw.WriteLine("v " + fornix.outSide.getUpver(i).X().ToString() + " "
                        + fornix.outSide.getUpver(i).Y().ToString() + " "
                        + fornix.outSide.getUpver(i).Z().ToString() + " ");
                }
                //生成底面外顶点
                for (i = 0; i < fornix.outSide.lowvers.Count; ++i)
                {
                    sw.WriteLine("v " + fornix.outSide.getDownver(i).X().ToString() + " "
                        + fornix.outSide.getDownver(i).Y().ToString() + " "
                        + fornix.outSide.getDownver(i).Z().ToString() + " ");
                }
            }
            sw.WriteLine("\n#内部点");
            foreach (Fornix fornix in fornixs)
            {
                //生成顶面内部点
                if (null == fornix.upFace.InnerPoints)
                    break;
                for (i = 0; i < fornix.upFace.InnerPoints.Count; ++i)
                {
                    fornix.upFace.InnerPoints.getVer(i).ID += (verNum - 1);
                    sw.WriteLine("v " + fornix.upFace.InnerPoints.getVer(i).X().ToString() + " "
                        + fornix.upFace.InnerPoints.getVer(i).Y().ToString() + " "
                        + fornix.upFace.InnerPoints.getVer(i).Z().ToString() + " ");
                }
            }

            foreach (Fornix fornix in fornixs)
            {
                //生成外侧面模型
                sw.WriteLine("\no " + fornix.name + "OUTSIDE\ng " + fornix.name + "OUTSIDE");
                for (i = 0; i < fornix.outSide.tris.Count; ++i)
                {
                    sw.WriteLine("f " + fornix.outSide.tris.Values[i].points.getVer(0).ID.ToString() + " "
                        + fornix.outSide.tris.Values[i].points.getVer(1).ID.ToString() + " "
                        + fornix.outSide.tris.Values[i].points.getVer(2).ID.ToString() + " ");
                }
                //生成顶面模型
                sw.WriteLine("\no " + fornix.name + "UPFACE\ng " + fornix.name + "UPFACE");
                for (i = 0; i < fornix.upFace.tris.Count; ++i)
                {
                    sw.WriteLine("f " + (fornix.upFace.tris.Values[i].points.getVer(0).ID).ToString() + " "
                        + (fornix.upFace.tris.Values[i].points.getVer(1).ID).ToString() + " "
                        + (fornix.upFace.tris.Values[i].points.getVer(2).ID).ToString() + " ");
                }
                //生成底面模型
                sw.WriteLine("\no " + fornix.name + "DOWNFACE\ng " + fornix.name + "DOMNFACE");
                for (i = 0; i < fornix.downFace.tris.Count; ++i)
                {
                    sw.WriteLine("f " + fornix.downFace.tris.Values[i].points.getVer(0).ID.ToString() + " "
                        + fornix.downFace.tris.Values[i].points.getVer(1).ID.ToString() + " "
                        + fornix.downFace.tris.Values[i].points.getVer(2).ID.ToString() + " ");
                }
            }
            sw.Close();
            sw.Dispose();
        }
    }
}
