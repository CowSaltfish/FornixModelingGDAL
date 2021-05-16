using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace FornixModelingGDAL
{
    class RasterWriter
    {
        public static void writeFornixObj(double[] datatable, int srcWidth, int srcHeight, string path)
        {
            StreamWriter sw = File.CreateText(path);
            int i, j;
            for (i = 0; i < srcWidth; ++i)
            {
                for (j = 0; j < srcHeight; ++j)
                    sw.Write(datatable[i * srcWidth + j].ToString() + " ");
                sw.Write("\n\n");
            }    
            sw.Close();
            sw.Dispose();
        }
    }
}
