@echo off
rem suppress warnings about using windows style paths.
set CYGWIN=nodosfilewarning
rem In case epics is built dynamic (i.e. dlls present in bin)
if exist "C:\Epics\extensions\bin\%EPICS_HOST_ARCH%\procServ.exe" (
	procServ --allow -n "NIKONKS-SYNC" -p pid.txt -L log.txt --logstamp -e "C:\WinPython-64bit-2.7.10.3\python-2.7.10.amd64\pythonw.exe" -i ^D^C 2005 nikonKsSync.py
) else (
	pythonw.exe nikonKsSync.py
)
