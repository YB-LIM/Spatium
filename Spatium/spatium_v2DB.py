from abaqusConstants import *
from abaqusGui import *
from kernelAccess import mdb, session
import os

thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)


###########################################################################
# Class definition
###########################################################################

class Spatium_v2DB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #

        AFXDataDialog.__init__(self, form, 'Spatium',
            self.OK|self.APPLY|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
            

        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
            

        applyBtn = self.getActionButton(self.ID_CLICKED_APPLY)
        applyBtn.setText('Apply')
            
        TabBook_1 = FXTabBook(p=self, tgt=None, sel=0,
            opts=TABBOOK_NORMAL,
            x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
            pt=DEFAULT_SPACING, pb=DEFAULT_SPACING)
        tabItem = FXTabItem(p=TabBook_1, text='Model setup', ic=None, opts=TAB_TOP_NORMAL,
            x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_1 = FXVerticalFrame(p=TabBook_1,
            opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
            x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
            pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        HFrame_1 = FXHorizontalFrame(p=TabItem_1, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        GroupBox_1 = FXGroupBox(p=HFrame_1, text='Geometry', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        ComboBox_1 = AFXComboBox(p=GroupBox_1, ncols=0, nvis=1, text='Type:                                               ', tgt=form.Is_PorousKw, sel=0)
        ComboBox_1.setMaxVisible(10)
        ComboBox_1.appendItem(text='Composite')
        ComboBox_1.appendItem(text='Porous')
        AFXTextField(p=GroupBox_1, ncols=12, labelText='Periodic Cell Size:                              ', tgt=form.LKw, sel=0)
        AFXTextField(p=GroupBox_1, ncols=12, labelText='Global Mesh Size:                              ', tgt=form.L_meshKw, sel=0)
        AFXTextField(p=GroupBox_1, ncols=12, labelText='Average Radius of Particle:                  ', tgt=form.r_avgKw, sel=0)
        AFXTextField(p=GroupBox_1, ncols=12, labelText='Standard Deviation of Particle Radius:   ', tgt=form.r_stdKw, sel=0)
        AFXTextField(p=GroupBox_1, ncols=12, labelText='Target Volume Fraction of Particle:       ', tgt=form.VoF_tarKw, sel=0)
        VFrame_7 = FXVerticalFrame(p=HFrame_1, opts=LAYOUT_FILL_X|LAYOUT_FILL_Y, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        GroupBox_5 = FXGroupBox(p=VFrame_7, text='Numerical Parameter', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        AFXTextField(p=GroupBox_5, ncols=12, labelText='Min. Distance Between Particle:          ', tgt=form.min_distanceKw, sel=0)
        AFXTextField(p=GroupBox_5, ncols=12, labelText='Max. Iteration for Particle Generation:  ', tgt=form.max_iterationsKw, sel=0)
        GroupBox_4 = FXGroupBox(p=VFrame_7, text='Boundary Condition', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        ComboBox_2 = AFXComboBox(p=GroupBox_4, ncols=0, nvis=1, text='Deformation Mode:    ', tgt=form.ModeKw, sel=0)
        ComboBox_2.setMaxVisible(10)
        ComboBox_2.appendItem(text='Uniaxial Tension')
        ComboBox_2.appendItem(text='Biaxial Tension')
        ComboBox_2.appendItem(text='Confined Compression')
        ComboBox_2.appendItem(text='Simple Shear')
        AFXTextField(p=GroupBox_4, ncols=12, labelText='Displacement:           ', tgt=form.DispKw, sel=0)
        fileName = os.path.join(thisDir, 'Main.png')
        icon = afxCreatePNGIcon(fileName)
        FXLabel(p=TabItem_1, text='', ic=icon)
        l = FXLabel(p=TabItem_1, text='Developed by Youngbin LIM', opts=JUSTIFY_LEFT)
        tabItem = FXTabItem(p=TabBook_1, text='Homogenization', ic=None, opts=TAB_TOP_NORMAL,
            x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_2 = FXVerticalFrame(p=TabBook_1,
            opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
            x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
            pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        VFrame_8 = FXVerticalFrame(p=TabItem_2, opts=LAYOUT_FILL_X|LAYOUT_FILL_Y, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        fileHandler = Spatium_v2DBFileHandler(form, 'Odb_Path', 'All files (*)')
        fileTextHf = FXHorizontalFrame(p=VFrame_8, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        # Note: Set the selector to indicate that this widget should not be
        #       colored differently from its parent when the 'Color layout managers'
        #       button is checked in the RSG Dialog Builder dialog.
        fileTextHf.setSelector(99)
        AFXTextField(p=fileTextHf, ncols=30, labelText='Odb file:                   ', tgt=form.Odb_PathKw, sel=0,
            opts=AFXTEXTFIELD_STRING|LAYOUT_CENTER_Y)
        icon = afxGetIcon('fileOpen', AFX_ICON_SMALL )
        FXButton(p=fileTextHf, text='	Select File\nFrom Dialog', ic=icon, tgt=fileHandler, sel=AFXMode.ID_ACTIVATE,
            opts=BUTTON_NORMAL|LAYOUT_CENTER_Y, x=0, y=0, w=0, h=0, pl=1, pr=1, pt=1, pb=1)
        fileHandler = Spatium_v2DBFileHandler(form, 'Output_Path', 'All files (*)')
        fileTextHf = FXHorizontalFrame(p=VFrame_8, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        # Note: Set the selector to indicate that this widget should not be
        #       colored differently from its parent when the 'Color layout managers'
        #       button is checked in the RSG Dialog Builder dialog.
        fileTextHf.setSelector(99)
        AFXTextField(p=fileTextHf, ncols=30, labelText='Output file:               ', tgt=form.Output_PathKw, sel=0,
            opts=AFXTEXTFIELD_STRING|LAYOUT_CENTER_Y)
        icon = afxGetIcon('fileOpen', AFX_ICON_SMALL )
        FXButton(p=fileTextHf, text='	Select File\nFrom Dialog', ic=icon, tgt=fileHandler, sel=AFXMode.ID_ACTIVATE,
            opts=BUTTON_NORMAL|LAYOUT_CENTER_Y, x=0, y=0, w=0, h=0, pl=1, pr=1, pt=1, pb=1)
        AFXTextField(p=VFrame_8, ncols=12, labelText='Nominal strain:          ', tgt=form.EngStrainKw, sel=0)
        ComboBox_3 = AFXComboBox(p=VFrame_8, ncols=0, nvis=1, text='Stress component:     ', tgt=form.S_CompKw, sel=0)
        ComboBox_3.setMaxVisible(10)
        ComboBox_3.appendItem(text='S11')
        ComboBox_3.appendItem(text='S22')
        ComboBox_3.appendItem(text='S33')
        ComboBox_3.appendItem(text='S12')
        ComboBox_3.appendItem(text='S13')
        ComboBox_3.appendItem(text='S23')
        FXCheckButton(p=VFrame_8, text='Plot curve on OK', tgt=form.Plot_FlagKw, sel=0)
        l2 = FXLabel(p=TabItem_2, text='*Note: Volume averaging is performed when both the ODB and output path are defined', opts=JUSTIFY_LEFT)        
        fileName = os.path.join(thisDir, 'Main2.png')
        icon = afxCreatePNGIcon(fileName)
        FXLabel(p=VFrame_8, text='', ic=icon)


###########################################################################
# Class definition
###########################################################################

class Spatium_v2DBFileHandler(FXObject):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form, keyword, patterns='*'):

        self.form = form
        self.patterns = patterns
        self.patternTgt = AFXIntTarget(0)
        exec('self.fileNameKw = form.%sKw' % keyword)
        self.readOnlyKw = AFXBoolKeyword(None, 'readOnly', AFXBoolKeyword.TRUE_FALSE)
        FXObject.__init__(self)
        FXMAPFUNC(self, SEL_COMMAND, AFXMode.ID_ACTIVATE, Spatium_v2DBFileHandler.activate)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def activate(self, sender, sel, ptr):

       fileDb = AFXFileSelectorDialog(getAFXApp().getAFXMainWindow(), 'Select a File',
           self.fileNameKw, self.readOnlyKw,
           AFXSELECTFILE_ANY, self.patterns, self.patternTgt)
       fileDb.setReadOnlyPatterns('*.odb')
       fileDb.create()
       fileDb.showModal()
