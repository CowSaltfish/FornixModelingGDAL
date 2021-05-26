//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;

//namespace FornixModelingGDAL
//{
//    class TDOverlayerAnalysisTriangularmesh : public TDconversionBase
//    {
//        public bool triangularmesh(LayerPathInfo originalPathInfo, LayerPathInfo targetPathInfo);
//        private bool readOriginalInfo(LayerPathInfo originalPathInfo, LayerPathInfo targetPathInfo, TranForm_enum enum_obj = tranform_save);
//        private void tranformToPoint(TranForm_enum enum_obj, OGRFeature poFeature) override;
//        private vector<OGRFeature> m_orgdatainfo;
//        private void calarminfo();
//        private void calarmVoronoiDiagram();

//        void TDOverlayerAnalysisTriangularmesh::calarmVoronoiDiagram()
//        {
//            VoronoiDiagramBuilder voronoiDiagramBuilder;
//            CoordinateArraySequence coords;
//            Envelope clipEnvelpoe;
//            std::for_each(m_orgdatainfo.begin(), m_orgdatainfo.end(), [&](OGRFeature *p) {
//                auto ppo = (OGRPoint *)(p->GetGeometryRef());
//                Coordinate coord(ppo->getX(), ppo->getY());
//                coords.add(coord);
//                clipEnvelpoe.expandToInclude(coord);
//            });
//            voronoiDiagramBuilder.setSites(coords);
//            voronoiDiagramBuilder.setClipEnvelope(&clipEnvelpoe);
//            auto vecgeom = voronoiDiagramBuilder.getDiagram(*GeometryFactory::getDefaultInstance());
//            for (int i = 0; i < vecgeom->getNumGeometries(); i++)
//            {
//                auto singlegeom = vecgeom->getGeometryN(i);
//                WKBWriter wkbwriter;
//                std::ostringstream ostr;
//                wkbwriter.write(*singlegeom, ostr);
//                OGRPolygon *armtemp = new OGRPolygon();
//                std::string  gstr = ostr.str();
//                int atoput;
//                armtemp->importFromWkb((unsigned char *)gstr.c_str(),-1, wkbVariantOldOgc, atoput);
//                auto armfeature = m_orgdatainfo[i]->Clone();
//                armfeature->SetGeometry(armtemp);
//                m_cachevector.push_back(armfeature);
//            }
//        }


//        void TDOverlayerAnalysisTriangularmesh::calarminfo()
//        {
//            double *dafx= new double[m_orgdatainfo.size()];
//            double *dafy = new double[m_orgdatainfo.size()];
//            int num = 0;
//            std::for_each(m_orgdatainfo.begin(), m_orgdatainfo.end(),[&](OGRFeature *p){
//            auto ppo = (OGRPoint *)(p->GetGeometryRef());
//            dafx[num] = ppo->getX();
//            dafy[num] = ppo->getY();
//            num++;
//            });
//            auto triangularmesh = GDALTriangulationCreateDelaunay(m_orgdatainfo.size(),dafx, dafy);
//            for (int i = 0; i < triangularmesh->nFacets; i++)
//            {
//                int indexone = triangularmesh->pasFacets[i].anVertexIdx[0];
//                int indextwo = triangularmesh->pasFacets[i].anVertexIdx[1];
//                int indexthress = triangularmesh->pasFacets[i].anVertexIdx[2];
//                auto point1 = m_orgdatainfo[indexone]->GetGeometryRef()->toPoint();
//                auto point2 = m_orgdatainfo[indextwo]->GetGeometryRef()->toPoint();
//                auto point3 = m_orgdatainfo[indexthress]->GetGeometryRef()->toPoint();
//                OGRTriangle *geom = new OGRTriangle(*point1, *point2, *point3);
//                auto armfeature = m_orgdatainfo[indexone]->Clone();
//                armfeature->SetGeometry(geom);
//                m_cachevector.push_back(armfeature);
//            }
//            delete[] dafx;
//            delete[] dafy;
//            GDALTriangulationFree(triangularmesh);
//        }

//        bool TDOverlayerAnalysisTriangularmesh::readOriginalInfo(LayerPathInfo originalPathInfo, LayerPathInfo targetPathInfo, TranForm_enum enum_obj)
//        {
//            (void)openOriginalInfo(originalPathInfo);
//            IF_NULL_RETRUN_FALSE(m_orgLayer)
//                OGRFeature *poFeature;
//            m_orgLayer->ResetReading();
//            OGRFeatureDefn * ogrfeaturedefn = m_orgLayer->GetLayerDefn();
//            pmoSourceSRS = m_orgLayer->GetSpatialRef();
//            while ((poFeature = m_orgLayer->GetNextFeature()) != NULL)
//            {
//                tranformToPoint(enum_obj, poFeature);
//                //OGRFeature::DestroyFeature(poFeature);
//            }
//            //calarminfo();
//            calarmVoronoiDiagram();
//            wiriteInfo(targetPathInfo, ogrfeaturedefn);
//            m_cachevector.clear();
//            if (NULL != m_pwriteDataset)
//            {
//                GDALClose(m_pwriteDataset);
//            }
//            if (NULL != m_preadDataset)
//            {
//                GDALClose(m_preadDataset);
//            }
//            return true;
//        }

//        bool TDOverlayerAnalysisTriangularmesh::triangularmesh(LayerPathInfo originalPathInfo, LayerPathInfo targetPathInfo)
//        {
//            return readOriginalInfo(originalPathInfo, targetPathInfo, transform_delaunaytriangulation);
//        }

//        void TDOverlayerAnalysisTriangularmesh::tranformToPoint(TranForm_enum enum_obj, OGRFeature *poFeature)
//        {

//            IF_NULL_RETRUN(poFeature)
//                auto pgeom = poFeature->GetGeometryRef();
//            if (NULL != pgeom && pgeom->getGeometryType() == wkbPoint)
//            {
//                m_orgdatainfo.push_back(poFeature);
//            }
//        }//`这里写代码片`
//    }
//}
