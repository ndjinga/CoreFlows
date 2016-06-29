from PyQt4 import QtGui, QtCore
from PyQt4.uic import loadUi
from utils import completeResPath

import CoreFlows as cfa

class MainCFWidget(QtGui.QTabWidget):
  def __init__(self):
    QtGui.QTabWidget.__init__(self)
    loadUi(completeResPath("MainCFWidget.ui"), self)
    self._python_dump = []
    
  def scanWidgets(self):
    print self.tabModel
    for k in self.__dict__:
      att = self.__dict__[k] 
      if isinstance(att, QtGui.QWidget):
        name = att.objectName()
        if name != "":
          if name.endswith("RadioButton"):
            assert(isinstance(att, QRadioButton))
            # parse name
            name=name[:len(name)-len("_RadioButton")]#On retire le suffixe _Radiobutton
            if name.startswith("Dim") :
		dictCF["spaceDim"]=name[3]

            if name==DriftModel or name==SinglePhase or name.endswith("Equation") or name.endswith("TwoFluid") :
		dictCF["ModelName"]=name
          elif name.endswith("doubleSpinBox"):
            assert(isinstance(att, QSpinBox))
            val = att.value()
            # parse name
            name=name[:len(name)-len("_doubleSpinBox")]#On retire le suffixe _doubleSpinBox
            dictCF[name]=val
          elif name.endswith('comboBox'):
            assert(isinstance(att, QComboBox))
            val = att.currentText()
            # parse name
            name=name[:len(name)-len("_comboBox")]#On retire le suffixe _comboBox
            if name==DM_1barOr155bar :
		pressureEstimate=val
            dictCF[name]=val 
    return dictCF  
          
  def onLaunchSimu(self):
    print "launching simulation"
    dictCF = self.scanWidgets()
    exec "myproblem = cfa.%s(cfa.%s,%s)" % (dictCF["ModelName"],dictCF["pressureEstimate"],dictCF["SpaceDim"])
    for k, v in dictCF.items():
        line = "cfa.set%s(%s)" % (k, v)
        #exec line
        self._python_dump.append(line)
    exec "myproblem.initialize()"
    exec "myproblem.run()"
    exec "myproblem.terminate()"
    
    # TODO use a helper object here.
