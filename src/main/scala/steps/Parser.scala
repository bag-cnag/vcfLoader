/*

case class rawTable(chrom : String, pos:Int, ID : String, ref :String ,alt : String, qual:String,filter:String,info : String, format:String,Sample : String)
    
def file_to_parquet(origin_path: String, destination : String, partition : String)=
{      //remove header
       val file = sc.textFile(origin_path).filter(line => !line.startsWith("#"))
       val raw_file = file.map(_.split("\t")).map(p => rawTable(p(0),p(1).trim.toInt,p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9))).toDF()
       raw_file.save(destination+"/sample="+partition)

}

file_to_parquet("/user/dpiscia/gvcf10bands/E000010.g.vcf.gz","/user/dpiscia/test/trio","E000010")
file_to_parquet("/user/dpiscia/gvcf10bands/E000036.g.vcf.gz","/user/dpiscia/test/trio","E000036")
file_to_parquet("/user/dpiscia/gvcf10bands/E000037.g.vcf.gz","/user/dpiscia/test/trio","E000037")


val sqlContext = new org.apache.spark.sql.hive.HiveContext(sc)

def reader(File : String)
//READ the FILE, name should pass as variable
   val file = sc.textFile("/user/admin/Chr22_VQSR110test.annot.snpEff.g.vcf.gz").filter(line => !line.startsWith("#"))
//val filtered_file = file.filter(line => line.split("\t")(1)=="51066632").take(1)
   // take sample columns out of main row
   //import org.apache.spark.sql.{SaveMode}
   //variants.save("variants",SaveMode.Append)
   for (i <- Range(9,119)) yield res.take(2)(1).split("\t")(i)
   
res.filter(x => !x.startsWith("#")).take(2).map(line => line.split("\t").zipWithIndex.filter(x => x._2 > 9 ).map(x => x._1.split(":")))


case class Sample(sample_id : String, gt : String, dp :Int, gq: Float)


case class Population(gp1_afr_af : Float,gp1_asn_af : Float, gp1_eur_af : Float, esp6500_aa: Float, esp6500_ea : Float, exac:Float, gmaf: Float, rd_freq : Float)

/*case class Effect(codon_change : String, amino_acid_change: String,amino_acid_length: String,effect : String, effect_impact : String ,exon_rank : String , functional_class : String,
                  gene_coding : String, gene_name : String, transcript_biotype : String, transcript_id: String,
                  gerp_plus_plus_rs : String,
                  cadd_phred : String, mt :String,mutationtaster_pred: String,phylop46way_placental: String, polyphen2_hvar_pred:String, polyphen2_hvar_score:String, 
                  sift_pred : String, sift_score : String, siphy_29way_pi :String )
*/

case class Effect(effect : String,effect_impact:String, functional_class : String,codon_change : String, amino_acid_change: String,amino_acid_length: String, gene_name : String,
    transcript_biotype : String,gene_coding : String, transcript_id: String,exon_rank : String, geno_type_number : Int)
case class Variant(chrom : String, pos: Int, ref : String, alt : String, indel : Boolean, rs : String, samples: List[Sample], population : List[Population], effects : List[Effect])
case class Variant(chrom : String, pos: Int, ref : String, alt : String, indel : Boolean, rs : String)


def parser(line : String): Variant=
{
   var fields = line.split("\t")
   var chrom = fields(0)
   var pos = fields(1).toInt
   var ref = fields(3)
   var alt= if (fields(4).contains("(") ) "multi"
     else fields(4)
   var indel= ref.length!=alt.length
   Variant(chrom,pos,ref,alt,indel,"")
   
}

def functional_parser(raw_line:String):Array[Effect]=
{
  val items=raw_line.split(",")
  items.map(item => {
    Effect(item.split("\\(")(0),
           item.split("\\(")(1),
           item.split("\\|")(1),
           item.split("\\|")(2),
           item.split("\\|")(3),
           item.split("\\|")(4),
           item.split("\\|")(5),
           item.split("\\|")(6),
           item.split("\\|")(7),
           item.split("\\|")(8),
           item.split("\\|")(9),
           item.split("\\|")(10).replace(")","").toInt)       
           
    
  })
}


val variants = file.map(line => parser(line)).toDF()

variants.registerTempTable("variants")
val some_variants = sqlContext.sql("SELECT * FROM variants limit 10")
some_variants.collect().foreach(println)
//variants.saveAsTable("variants_prova")

variants.saveAsParquetFile("variants_prova.parquet")

//for effects input

(7).split("EFF=")(1)

val effects = parser_functional(res(0).split("\t")(7).split("EFF=")(1))
def prediction_parser(info :String)={
  val null_option=""
  val result = info.split(";").map(_ split "=") collect { case Array(k, v) => (k, v) } toMap
  val sift_pred = result.getOrElse("SIFT_pred",null_option)
  val sift_score = result.getOrElse("SIFT_score",null_option)
  val Polyphen2_HVAR_pred = result.getOrElse("Polyphen2_HVAR_pred",null_option)
  val pp2 = result.getOrElse("pp2",null_option)
  val Polyphen2_HVAR_score = result.getOrElse("Polyphen2_HVAR_score",null_option)
  val MutationTaster_pred = result.getOrElse("MutationTaster_pred",null_option)
  val mt = result.getOrElse("mt",null_option)
  val phyloP46way_placental = result.getOrElse("phyloP46way_placental",null_option)
  val GERP = result.getOrElse("GERP",null_option)
  val SiPhy_29way_pi = result.getOrElse("SiPhy_29way_pi",null_option)
  val CADD_phred = result.getOrElse("CADD_phred",null_option)
  (sift_pred,sift_score,Polyphen2_HVAR_pred,pp2,Polyphen2_HVAR_score,MutationTaster_pred,mt,phyloP46way_placental,GERP,SiPhy_29way_pi,CADD_phred)
}


def global_parser(line :String)={
  //if eff is not present it will crash
  val fields=line.split("\t")
  val chrom = fields(0)
  val pos = fields(1).toInt
  val ref = fields(3)
  val alt=  fields(4).split(",")
  val indel= ref.length!=alt.length  
  val effects = parser_functional(fields(7).split("EFF=")(1))
  println("fields 2 is "+fields(2))
  val predictions = prediction_parser(fields(2))
  (chrom,pos,ref,alt,indel,effects,predictions)
   
}
val res = file.map(line => global_parser(line)).take(10)
//val filtered = file.filter(line => line.split("\t")(4).contains(",")).filter(line=> line.split("\t")(2).contains("SIFT_pred")).take(1)

*/