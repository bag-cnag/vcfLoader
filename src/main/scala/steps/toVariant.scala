package steps


        
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

  annotations.join(samples, annotations("pos") === samples("_1") && annotations("ref") === samples("_2") && annotations("alt") === samples("_3"), "right")
    .select("pos","ref","alt","rs","indel","_6","_c5","_c6","_c7")
    .withColumnRenamed("_6","samples")
    .withColumnRenamed("_c5","effs")
    .withColumnRenamed("_c6","populations")
    .withColumnRenamed("_c7","predictions")
    .save(destination+"/chrom="+chromList)//+"/band="+banda._2.toString)

}
}