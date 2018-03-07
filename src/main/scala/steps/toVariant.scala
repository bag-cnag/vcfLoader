package steps

import org.apache.spark.sql.functions.udf

    case class VariantModel(pos:Int,ref:String,alt:String,rs:String,indel:Boolean,
         samples: collection.mutable.WrappedArray[Map[String,String]],
        effs: collection.mutable.WrappedArray[Map[String,String]],
        populations:collection.mutable.WrappedArray[Map[String,String]],
        predictions:collection.mutable.WrappedArray[Map[String,String]])
        
        object toVariant {
def main(sc :org.apache.spark.SparkContext, Samples:org.apache.spark.sql.DataFrame, Annotations:org.apache.spark.sql.DataFrame,
        destination: String,
    chromList : String, 
    banda : (Int,Int))={
val sqlContext = new org.apache.spark.sql.SQLContext(sc)

// this is used to implicitly convert an RDD to a DataFrame.
import sqlContext.implicits._

//pos _2,ref_3,alt_4,rs_5,indel_6, smaples_7
val samples = Samples
    .where(Samples("chrom")===chromList.toInt)
    //.where(Samples("band") ===banda._2)
    val annotations = Annotations
      .where(Annotations("chrom")===chromList.toInt)

  val freq = udf((array: scala.collection.mutable.WrappedArray[Map[String,String]]) => {
    val temp = array.map(variant=> variant.getOrElse("gt","0/0"));
    temp.toList.flatMap(x=>x.split("/")).map(x=> x.toInt).sum/((2*temp.length).toFloat) })
  val freqCounts = udf((array: scala.collection.mutable.WrappedArray[Map[String,String]]) => {
    val temp = array.map(variant=> variant.getOrElse("gt","0/0"));
    temp.toList.length*2 })

  annotations.join(samples, annotations("pos2") === samples("pos") && annotations("ref2") === samples("ref") && annotations("alt2") === samples("alt"), "left")
    .select("pos","ref","alt","indel","samples","effs","populations","predictions")
    .withColumn("freqInt",freq(samples("samples")))
    .withColumn("freqCounts",freqCounts(samples("samples")))
    .write.parquet(destination+"/chrom="+chromList)//+"/band="+banda._2.toString)


}
}
