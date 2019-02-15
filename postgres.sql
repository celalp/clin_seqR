-- ---
-- Globals
-- ---

-- SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
-- SET FOREIGN_KEY_CHECKS=0;

-- ---
-- Table 'users'
-- 
-- ---

DROP TABLE IF EXISTS `users`;
		
CREATE TABLE `users` (
  `userid` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `username` VARCHAR NOT NULL DEFAULT 'NULL',
  `projectid` INTEGER NOT NULL DEFAULT NULL,
  PRIMARY KEY (`userid`)
);

-- ---
-- Table 'projects'
-- 
-- ---

DROP TABLE IF EXISTS `projects`;
		
CREATE TABLE `projects` (
  `projectid` INTEGER NOT NULL AUTO_INCREMENT DEFAULT NULL,
  `projectname` VARCHAR NULL DEFAULT NULL,
  `description` MEDIUMTEXT NULL DEFAULT NULL,
  PRIMARY KEY (`projectid`)
);

-- ---
-- Table 'samples'
-- 
-- ---

DROP TABLE IF EXISTS `samples`;
		
CREATE TABLE `samples` (
  `sampleid` INTEGER NOT NULL AUTO_INCREMENT DEFAULT NULL,
  `samplename` VARCHAR NULL DEFAULT NULL,
  `projectid` INTEGER NULL DEFAULT NULL,
  `tissue` VARCHAR NOT NULL DEFAULT 'NULL',
  `description` MEDIUMTEXT NULL DEFAULT NULL,
  `files` VARCHAR NULL DEFAULT NULL,
  PRIMARY KEY (`sampleid`),
KEY (`sampleid`)
);

-- ---
-- Table 'genomes'
-- 
-- ---

DROP TABLE IF EXISTS `genomes`;
		
CREATE TABLE `genomes` (
  `genomeid` INTEGER NOT NULL AUTO_INCREMENT DEFAULT NULL,
  `genomename` VARCHAR NOT NULL DEFAULT 'NULL',
  `organism` VARCHAR NULL DEFAULT NULL,
  PRIMARY KEY (`genomeid`)
);

-- ---
-- Table 'data'
-- 
-- ---

DROP TABLE IF EXISTS `data`;
		
CREATE TABLE `data` (
  `dataid` INTEGER NOT NULL AUTO_INCREMENT DEFAULT NULL,
  `name` VARCHAR NULL DEFAULT NULL,
  `source` VARCHAR NULL DEFAULT NULL,
  `version` VARCHAR NULL DEFAULT NULL,
  `genomeid` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`dataid`)
);

-- ---
-- Table 'tools'
-- 
-- ---

DROP TABLE IF EXISTS `tools`;
		
CREATE TABLE `tools` (
  `toolid` INTEGER NOT NULL AUTO_INCREMENT DEFAULT NULL,
  `name` VARCHAR NULL DEFAULT NULL,
  `version` VARCHAR NULL DEFAULT NULL,
  PRIMARY KEY (`toolid`)
);

-- ---
-- Table 'runs'
-- 
-- ---

DROP TABLE IF EXISTS `runs`;
		
CREATE TABLE `runs` (
  `runid` INTEGER NOT NULL AUTO_INCREMENT DEFAULT NULL,
  `pipelineid` INTEGER NULL DEFAULT NULL,
  `sampleid` INTEGER NULL DEFAULT NULL,
  `date` INTEGER NULL DEFAULT NULL,
  `genomeid` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`runid`)
);

-- ---
-- Table 'genes'
-- 
-- ---

DROP TABLE IF EXISTS `genes`;
		
CREATE TABLE `genes` (
  `geneid` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `genename` VARCHAR NULL DEFAULT NULL,
  `genesymbol` VARCHAR NULL DEFAULT NULL,
  `genomeid` INTEGER NULL DEFAULT NULL,
  `start` INTEGER NULL DEFAULT NULL,
  `end` INTEGER NULL DEFAULT NULL,
  `chr` INTEGER NULL DEFAULT NULL,
  `strand` INTEGER NULL DEFAULT NULL,
  `info` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`geneid`)
);

-- ---
-- Table 'transcripts'
-- 
-- ---

DROP TABLE IF EXISTS `transcripts`;
		
CREATE TABLE `transcripts` (
  `txid` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `txname` VARCHAR NULL DEFAULT NULL,
  `geneid` INTEGER NULL DEFAULT NULL,
  `start` INTEGER NULL DEFAULT NULL,
  `end` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`txid`)
);

-- ---
-- Table 'exons'
-- 
-- ---

DROP TABLE IF EXISTS `exons`;
		
CREATE TABLE `exons` (
  `exonid` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `exnname` VARCHAR NULL DEFAULT NULL,
  `txid` INTEGER NULL DEFAULT NULL,
  `start` INTEGER NULL DEFAULT NULL,
  `end` INTEGER NULL DEFAULT NULL,
  `rank` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`exonid`)
);

-- ---
-- Table 'junctions'
-- 
-- ---

DROP TABLE IF EXISTS `junctions`;
		
CREATE TABLE `junctions` (
  `juncid` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `txid` INTEGER NULL DEFAULT NULL,
  `ex1id` INTEGER NULL DEFAULT NULL,
  `ex2id` INTEGER NULL DEFAULT NULL,
  `start` INTEGER NULL DEFAULT NULL,
  `end` INTEGER NULL DEFAULT NULL,
  `motif` CHAR NULL DEFAULT NULL,
  `annotated` BINARY NULL DEFAULT NULL,
  PRIMARY KEY (`juncid`)
);

-- ---
-- Table 'variants'
-- 
-- ---

DROP TABLE IF EXISTS `variants`;
		
CREATE TABLE `variants` (
  `variantid` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `chr` INTEGER NULL DEFAULT NULL,
  `start` INTEGER NULL DEFAULT NULL,
  `ref` VARCHAR NULL DEFAULT NULL,
  `alt` VARCHAR NULL DEFAULT NULL,
  `genomeid` INTEGER NULL DEFAULT NULL,
  `info` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`variantid`)
);

-- ---
-- Table 'sample_variants'
-- 
-- ---

DROP TABLE IF EXISTS `sample_variants`;
		
CREATE TABLE `sample_variants` (
  `id` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `runid` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`id`)
);

-- ---
-- Foreign Keys 
-- ---

ALTER TABLE `users` ADD FOREIGN KEY (projectid) REFERENCES `projects` (`projectid`);
ALTER TABLE `samples` ADD FOREIGN KEY (projectid) REFERENCES `projects` (`projectid`);
ALTER TABLE `data` ADD FOREIGN KEY (genomeid) REFERENCES `genomes` (`genomeid`);
ALTER TABLE `runs` ADD FOREIGN KEY (sampleid) REFERENCES `samples` (`sampleid`);
ALTER TABLE `runs` ADD FOREIGN KEY (genomeid) REFERENCES `genomes` (`genomeid`);
ALTER TABLE `genes` ADD FOREIGN KEY (genomeid) REFERENCES `genomes` (`genomeid`);
ALTER TABLE `transcripts` ADD FOREIGN KEY (geneid) REFERENCES `genes` (`geneid`);
ALTER TABLE `exons` ADD FOREIGN KEY (txid) REFERENCES `transcripts` (`txid`);
ALTER TABLE `junctions` ADD FOREIGN KEY (txid) REFERENCES `transcripts` (`txid`);
ALTER TABLE `junctions` ADD FOREIGN KEY (ex1id) REFERENCES `exons` (`exonid`);
ALTER TABLE `junctions` ADD FOREIGN KEY (ex2id) REFERENCES `exons` (`exonid`);
ALTER TABLE `variants` ADD FOREIGN KEY (variantid) REFERENCES `sample_variants` (`id`);
ALTER TABLE `variants` ADD FOREIGN KEY (genomeid) REFERENCES `genomes` (`genomeid`);
ALTER TABLE `sample_variants` ADD FOREIGN KEY (runid) REFERENCES `runs` (`runid`);

-- ---
-- Table Properties
-- ---

-- ALTER TABLE `users` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `projects` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `samples` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `genomes` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `data` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `tools` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `runs` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `genes` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `transcripts` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `exons` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `junctions` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `variants` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;
-- ALTER TABLE `sample_variants` ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- ---
-- Test Data
-- ---

-- INSERT INTO `users` (`userid`,`username`,`projectid`) VALUES
-- ('','','');
-- INSERT INTO `projects` (`projectid`,`projectname`,`description`) VALUES
-- ('','','');
-- INSERT INTO `samples` (`sampleid`,`samplename`,`projectid`,`tissue`,`description`,`files`) VALUES
-- ('','','','','','');
-- INSERT INTO `genomes` (`genomeid`,`genomename`,`organism`) VALUES
-- ('','','');
-- INSERT INTO `data` (`dataid`,`name`,`source`,`version`,`genomeid`) VALUES
-- ('','','','','');
-- INSERT INTO `tools` (`toolid`,`name`,`version`) VALUES
-- ('','','');
-- INSERT INTO `runs` (`runid`,`pipelineid`,`sampleid`,`date`,`genomeid`) VALUES
-- ('','','','','');
-- INSERT INTO `genes` (`geneid`,`genename`,`genesymbol`,`genomeid`,`start`,`end`,`chr`,`strand`,`info`) VALUES
-- ('','','','','','','','','');
-- INSERT INTO `transcripts` (`txid`,`txname`,`geneid`,`start`,`end`) VALUES
-- ('','','','','');
-- INSERT INTO `exons` (`exonid`,`exnname`,`txid`,`start`,`end`,`rank`) VALUES
-- ('','','','','','');
-- INSERT INTO `junctions` (`juncid`,`txid`,`ex1id`,`ex2id`,`start`,`end`,`motif`,`annotated`) VALUES
-- ('','','','','','','','');
-- INSERT INTO `variants` (`variantid`,`chr`,`start`,`ref`,`alt`,`genomeid`,`info`) VALUES
-- ('','','','','','','');
-- INSERT INTO `sample_variants` (`id`,`runid`) VALUES
-- ('','');