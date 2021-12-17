workflow smCounterWF {

  String name
  String captureID
  File inputFasta1
  File inputFasta2
  File paramsFile
  File primerFile
  File roiFile
  File qiaseqCont
  File qiaseqRun

  call smCounter {
    input:
      SampleName=name,
      CaptureID=captureID,
      InputFasta1=inputFasta1,
      InputFasta2=inputFasta2,
      ParamsFile=paramsFile,
      PrimerFile=primerFile,
      RoiFile=roiFile,
      QiaseqCont=qiaseqCont,
      QiaseqRun=qiaseqRun
  }
  
}

task smCounter {

  String SampleName
  String CaptureID
  File InputFasta1
  File InputFasta2
  File ParamsFile
  File PrimerFile
  File RoiFile
  File QiaseqCont
  File QiaseqRun

  command {
    cat ${ParamsFile} | \
    sed "s/insert_sample_name/${SampleName}/g" | \
    sed "s/insert_capture_id/${CaptureID}/g" \
    > run_sm_counter_v2.params.txt
  }
  output {

  }
  
}
