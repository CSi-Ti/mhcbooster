
@echo off
set PYTHON="C:\Users\rw762\AppData\Local\miniconda3\python"
set JAVA="D:\software\fragpipe\jre\bin\java.exe"
set MSFragger="D:\software\fragpipe\tools\MSFragger-4.1\MSFragger-4.1.jar"
set MSFragger_SPLIT="D:\software\fragpipe\tools\msfragger_pep_split.py"
set N_SPLIT=4
set PARAMS="E:\test_calmzml\fragger.params"

@REM %JAVA% -jar -Dfile.encoding=UTF-8 -Xmx30G %MSFragger% %PARAMS% E:\test\JY_Class1_1M_DDA_60min_Slot1-10_1_541.d

%PYTHON% %MSFragger_SPLIT% %N_SPLIT% "%JAVA% -jar -Dfile.encoding=UTF-8 -Xmx30G" %MSFragger% %PARAMS% E:\test_calmzml\JY_Class1_1M_DDA_60min_Slot1-10_1_541.d E:\test_calmzml\JY_Class1_10M_DDA_60min_Slot1-11_1_545.d E:\test_calmzml\JY_Class1_25M_DDA_60min_Slot1-12_1_552.d
