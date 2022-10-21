# NOTE: This file will get overwritten when updating Qt Creator.
#
# To add dumpers that don't get overwritten, copy this file here
# to a safe location outside the Qt Creator installation and
# make this location known to Qt Creator using the Debugger >
# Locals & Expressions > Extra Debugging Helpers setting.

## original path:
##   c:\Qt\Tools\QtCreator\share\qtcreator\debugger\personaltypes.py
## https://doc.qt.io/qtcreator/creator-debugging-helpers.html#adding-custom-debugging-helpers
## Edit >  Preferences / Debugging Tool > Locals & Expressions > Extra Debugging Helpers:
##   C:\home\Development\Mandelbrot\personaltypes.py
## Debugging Helper Customization:
##   python theDumper.reloadDumpers({})


from dumper import Children, SubItem, UnnamedSubItem, DumperBase
from utils import DisplayFormat, TypeCode

######################## Your code below #######################

#so, I guess
#putValue, putType is for raw data
#putItem is for value[...] data

def qdump__MandelMath__worker_multi_double(d, value):
    d.putExpandable()
    if d.isExpanded():
      with Children(d):
        with SubItem(d, "(_ntype)"):
          d.putItem(value["_ntype"])
        #with SubItem(d, "storage_"):
        #  d.putItem(value["storage_"])
        with SubItem(d, "allocator"):
          d.putItem(value["allocator"])
        with SubItem(d, "capac"):
          #d.putValue(value["allocator"]["capacity"])
          #d.putItem(value["allocator"]["capacity"])
          d.putValue(value["allocator"]["capacity"].integer())
        with SubItem(d, "storage_"):
          #d.putValue(value["storage_"])
          #d.putArrayData(value["storage_"].address(), 4, d.lookupType('double'))
          d.putArrayData(value["storage_"].address(), 4+value["allocator"]["capacity"].integer(), d.lookupType('double'))          
        with SubItem(d, "test-d1"):
          d.putValue("Test-d")
    ##d.putArrayData(value["storage_"].address(), value["allocator"]["capacity"], d.lookupType('double'))
    #d.putExpandable()
    #if d.isExpanded():
    #d.putArrayData(value["storage_"].pointer(), 4, d.lookupType('double'))
    #d.putValue(value["storage_"].dereference())
    #d.putExpandable()
    #if d.isExpanded():
    #  with Children(d):
    #    with SubItem(d, "kuk"):
    #      d.putItem(value["_ntype"])
      #  with SubItem(d, "1"):
      #    d.putItem(value["storage_"][1])
      #d.putArrayData(value["storage_"].address(), 4, d.lookupType('double'))

def qdump__MandelMath__coXmplex(d, value):
    worker_type=value["allocator"]["worker"].type.name
    if worker_type=="MandelMath::worker_multi_double*":
      re_qt=value["re"]["asf64"].dereference()
      re_str="%.4g" % (re_qt.nativeValue, )   #value["re"]["asf64"].dereference().split("d")
      im_qt=value["im"]["asf64"].dereference()
      im_str="%+.4g" % (im_qt.nativeValue, )   #value["im"]["asf64"].dereference().split("d")
    elif worker_type=="MandelMath::worker_multi_float128*":
      re_qt=value["re"]["asf128"].dereference()
      re_str="%.4g" % (re_qt.nativeValue, )
      im_qt=value["im"]["asf128"].dereference()
      im_str="%+.4g" % (im_qt.nativeValue, )
    elif worker_type=="MandelMath::worker_multi_ddouble*":
      re_qt=value["re"]["asdd"].dereference()["hi"]
      re_str="%.4g" % (re_qt.nativeValue, )
      im_qt=value["im"]["asdd"].dereference()["hi"]
      im_str="%+.4g" % (im_qt.nativeValue, )
    elif worker_type=="MandelMath::worker_multi_real642*":
      re_qt=value["re"]["as642"].dereference()["val1"]
      re_str="%.4g" % (re_qt.nativeValue, )
      im_qt=value["im"]["as642"].dereference()["val1"]
      im_str="%+.4g" % (im_qt.nativeValue, )
    else:
      re_qt=None
      re_str="[!"+worker_type+"]"
      im_qt=None
      im_str="[?]"

    d.putExpandable()
    if d.isExpanded():
      with Children(d):
        with SubItem(d, "allocator"):
          d.putItem(value["allocator"])
        #d.putValue(str(value["allocator"]["worker"]))
        #with SubItem(d, "worker"):
        #  d.putValue(worker_type)#value["allocator"]["worker"])
        #  d.putType("def")#worker_type)
        with SubItem(d, "re"):
          d.putItem(value["re"]) #I'm probably hacking something, order matters here
          if re_qt is not None:
            d.putItem(re_qt)
        with SubItem(d, "im"):
          d.putItem(value["im"])
          if im_qt is not None:
            d.putItem(im_qt)
    d.putValue("%s %si" % (re_str, im_str, ))

def qdump__MandelMath__nuXmber(d, value):
    worker_type=value["allocator"]["worker"].type.name
    if worker_type=="MandelMath::worker_multi_double*":
      val_qt=value["ptr"]["asf64"].dereference()
      val_str="%.8g" % (val_qt.nativeValue, )
    elif worker_type=="MandelMath::worker_multi_float128*":
      val_qt=value["ptr"]["asf128"].dereference()
      val_str="%.8g" % (val_qt.nativeValue, )
    elif worker_type=="MandelMath::worker_multi_ddouble*":
      val_qt=value["ptr"]["asdd"].dereference()["hi"]
      val_str="%.8g" % (val_qt.nativeValue, )
    elif worker_type=="MandelMath::worker_multi_real642*":
      val_qt=value["ptr"]["as642"].dereference()["val1"]
      val_str="%.8g" % (val_qt.nativeValue, )
    else:
      val_qt=None
      val_str="[!"+worker_type+"]"

    d.putExpandable()
    if d.isExpanded():
      with Children(d):
        with SubItem(d, "allocator"):
          d.putItem(value["allocator"])
        with SubItem(d, "ptr"):
          d.putItem(value["ptr"]) #I'm probably hacking something, order matters here
          if val_qt is not None:
            d.putItem(val_qt)
    d.putValue(val_str)
