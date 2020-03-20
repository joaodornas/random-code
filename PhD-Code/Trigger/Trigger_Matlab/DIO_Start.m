out=daqhwinfo('agilentu2300');
dio=digitalio('agilentu2300',0)
dio

addline(dio,0:15,'in')
set(dio,'TimerFcn',@daqcallback);
set(dio,'TimerPeriod',5.0);
start(dio)
dio
stop(dio)

delete(daqfind);
delete(dio);
clear dio;