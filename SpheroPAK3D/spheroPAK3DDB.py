from abaqusConstants import *
from abaqusGui import *
from kernelAccess import mdb, session
import os

thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)


###########################################################################
# Class definition
###########################################################################

class SpheroPAK3DDB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #

        AFXDataDialog.__init__(self, form, 'SpheroPAK3D',
            self.OK|self.APPLY|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
            

        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
            

        applyBtn = self.getActionButton(self.ID_CLICKED_APPLY)
        applyBtn.setText('Apply')
            
        HFrame_6 = FXHorizontalFrame(p=self, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        GroupBox_6 = FXGroupBox(p=HFrame_6, text='Geometry', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_6 = AFXVerticalAligner(p=GroupBox_6, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        ComboBox_3 = AFXComboBox(p=VAligner_6, ncols=0, nvis=1, text='Type:', tgt=form.Is_PorousKw, sel=0)
        ComboBox_3.setMaxVisible(10)
        ComboBox_3.appendItem(text='Composite')
        ComboBox_3.appendItem(text='Porous')
        AFXTextField(p=VAligner_6, ncols=11, labelText='Periodic Cell Size:', tgt=form.LKw, sel=0)
        AFXTextField(p=VAligner_6, ncols=11, labelText='Global Mesh Size:', tgt=form.L_meshKw, sel=0)
        AFXTextField(p=VAligner_6, ncols=11, labelText='Average Radius of Particle:', tgt=form.r_avgKw, sel=0)
        AFXTextField(p=VAligner_6, ncols=11, labelText='Standard Deviation of Partical Radius:   ', tgt=form.r_stdKw, sel=0)
        AFXTextField(p=VAligner_6, ncols=11, labelText='Target Volume Fraction of Particle:', tgt=form.VoF_tarKw, sel=0)
        VFrame_6 = FXVerticalFrame(p=HFrame_6, opts=LAYOUT_FILL_X|LAYOUT_FILL_Y, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        GroupBox_17 = FXGroupBox(p=VFrame_6, text='Numerical Parameter', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_9 = AFXVerticalAligner(p=GroupBox_17, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=VAligner_9, ncols=11, labelText='Min. Distance Between Particle:', tgt=form.min_distanceKw, sel=0)
        AFXTextField(p=VAligner_9, ncols=11, labelText='Max. iteration for Particle Generation:  ', tgt=form.max_iterationsKw, sel=0)
        GroupBox_16 = FXGroupBox(p=VFrame_6, text='Boundary Condition', opts=FRAME_GROOVE|LAYOUT_FILL_Y)
        VAligner_10 = AFXVerticalAligner(p=GroupBox_16, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        ComboBox_5 = AFXComboBox(p=VAligner_10, ncols=0, nvis=1, text='Deformation Mode:    ', tgt=form.ModeKw, sel=0)
        ComboBox_5.setMaxVisible(10)
        ComboBox_5.appendItem(text='Uniaxial Tension')
        ComboBox_5.appendItem(text='Biaxial Tension')
        ComboBox_5.appendItem(text='Confined Compression')
        ComboBox_5.appendItem(text='Simple Shear')
        AFXTextField(p=VAligner_10, ncols=11, labelText='Displacement:   ', tgt=form.DispKw, sel=0)
        fileName = os.path.join(thisDir, 'Main.png')
        icon = afxCreatePNGIcon(fileName)
        FXLabel(p=self, text='', ic=icon)
        l = FXLabel(p=self, text='  Developed by Youngbin LIM', opts=JUSTIFY_LEFT)
