DROP TABLE IF EXISTS `pubmedlinks`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pubmedlinks` (
  `idLinks` int(11) NOT NULL AUTO_INCREMENT,
  `PMID` varchar(45) DEFAULT NULL,
  `PUBMED_Response` TEXT DEFAULT NULL,
  `PUBMED_Link_doi` TEXT DEFAULT NULL,
  `PUBMED_Link_FullText` TEXT DEFAULT NULL,
  `PUBMED_Link_PMC` TEXT DEFAULT NULL,
  `JOURNAL_Response_doi` TEXT DEFAULT NULL,
  `STATUS_doi` TEXT DEFAULT NULL,
  `S_vector_doi` varchar(45) DEFAULT NULL,
  `JOURNAL_Response_FullText` TEXT DEFAULT NULL,
  `STATUS_FullText` TEXT DEFAULT NULL,
  `S_vector_FullText` varchar(45) DEFAULT NULL,
  `JOURNAL_Response_PMC` TEXT DEFAULT NULL,
  `STATUS_PMC` TEXT DEFAULT NULL,
  `S_vector_PMC` varchar(45) DEFAULT NULL,
  `QueryLabelVersion` varchar(45) DEFAULT NULL,
  `PDFSTATUS` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`idLinks`)
) ENGINE=InnoDB AUTO_INCREMENT=0 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;