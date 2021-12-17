workflow UMI_dedup {

  File picard
  File inputBam
  File dedupMetrics
  String name

  call nameSort {
    input:
      SampleName=name,
      Picard=picard,
      InputBam=inputBam,
      DedupMetrics=dedupMetrics
  }
  call collapseUMI {
    input:
      SampleName=name,
      Picard=picard,
      DedupMetrics=dedupMetrics,
      NamesortedBam=nameSort.NamesortedBam
  }
  call coordSort {
    input:
      SampleName=name,
      Picard=picard,
      DedupBam=collapseUMI.DedupBam
  }
  call indexBam {
    input:
      SampleName=name,
      Picard=picard,
      CoordsortedBam=coordSort.CoordsortedBam,
      NamesortedBam=nameSort.NamesortedBam,
      DedupBam=collapseUMI.DedupBam
  }
  
}

task nameSort {

  File Picard
  File InputBam
  File DedupMetrics
  String SampleName

  command {
    samtools sort -n -o ${SampleName}.align.sorted.by.query.bam ${InputBam}
  }
  output {
    File NamesortedBam = "${SampleName}.align.sorted.by.query.bam"
  }
  
}

task collapseUMI {

  File Picard
  File DedupMetrics
  String SampleName
  File NamesortedBam

  command {
    java -jar ${Picard} MarkDuplicates \
      -I ${NamesortedBam} \
      -O ${SampleName}.dedup.bam \
      -M ${DedupMetrics} \
      -ASSUME_SORT_ORDER queryname \
      -BARCODE_TAG "mi" \
      -REMOVE_DUPLICATES True
  }
  output {
    File DedupBam = "${SampleName}.dedup.bam"
  }
  
}

task coordSort {

  File Picard
  File DedupBam
  String SampleName

  command {
    samtools sort -o ${SampleName}.dedup.sorted.by.coord.bam ${DedupBam}
  }
  output {
    File CoordsortedBam = "${SampleName}.dedup.sorted.by.coord.bam"
  }
  
}

task indexBam {

  File Picard
  File CoordsortedBam
  File NamesortedBam
  File DedupBam
  String Base = basename(CoordsortedBam)
  String SampleName

  command {
    mv ${CoordsortedBam} ${Base}
    samtools index ${Base}
    mv ${Base}.bai ${SampleName}.dedup.sorted.by.coord.bam.bai
    rm ${SampleName}.align.sorted.by.query.bam ${SampleName}.dedup.bam
  }
  output {
    File IndexedBam = "${SampleName}.dedup.sorted.by.coord.bam.bai"
  }
  
}
