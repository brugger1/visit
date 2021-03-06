// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <QvisLine2DInterface.h>

#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QLayout>
#include <QToolTip>

#include <AnnotationObject.h>

#include <QvisColorButton.h>
#include <QvisLineWidthWidget.h>
#include <QvisOpacitySlider.h>
#include <QvisScreenPositionEdit.h>

// ****************************************************************************
// Method: QvisLine2DInterface::QvisLine2DInterface
//
// Purpose: 
//   Constructor for the QvisLine2DInterface class.
//
// Arguments:
//   parent : This widget's parent widget.
//   name   : The name of this widget.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:31:16 PDT 2004
//
// Modifications:
//   Brad Whitlock, Tue Jun 28 13:36:41 PST 2005
//   Added code to make tool tips for the start, end coordinates.
//
//   Brad Whitlock, Tue Apr  8 09:27:26 PDT 2008
//   Support for internationalization.
//
//   Brad Whitlock, Mon Jul 21 10:39:10 PDT 2008
//   Qt 4.
//
//   Kathleen Biagas, Mon Jul 13 13:03:45 PDT 2015
//   Add useForegroundColor, colorlabel.
//
// ****************************************************************************

QvisLine2DInterface::QvisLine2DInterface(QWidget *parent) :
    QvisAnnotationObjectInterface(parent)
{
    // Set the title of the group box.
    this->setTitle(GetName());

    QGridLayout *cLayout = new QGridLayout(0);
    topLayout->addLayout(cLayout);
    cLayout->setSpacing(10);

    int row = 0;
    // Add controls for the start position
    positionStartEdit = new QvisScreenPositionEdit(this);
    connect(positionStartEdit, SIGNAL(screenPositionChanged(double, double)),
            this, SLOT(positionStartChanged(double, double)));
    QLabel *startLabel = new QLabel(tr("Start"), this);
    QString startTip(tr("Start of line in screen coordinates [0,1]"));
    startLabel->setToolTip(startTip);
    cLayout->addWidget(positionStartEdit, row, 1, 1, 3);
    cLayout->addWidget(startLabel, row, 0);
    ++row;

    // Add controls for the end position
    positionEndEdit = new QvisScreenPositionEdit(this);
    connect(positionEndEdit, SIGNAL(screenPositionChanged(double, double)),
            this, SLOT(positionEndChanged(double, double)));
    QLabel *endLabel = new QLabel(tr("End"), this);
    QString endTip(tr("End of line in screen coordinates [0,1]"));
    endLabel->setToolTip(endTip);
    cLayout->addWidget(positionEndEdit, row, 1, 1, 3);
    cLayout->addWidget(endLabel, row, 0);
    ++row;
   
    // Add controls for width.
    widthWidget = new QvisLineWidthWidget(0, this);
    connect(widthWidget, SIGNAL(lineWidthChanged(int)),
            this, SLOT(widthChanged(int)));
    cLayout->addWidget(widthWidget, row, 1);
    cLayout->addWidget(new QLabel(tr("Width"), this), row, 0);
    ++row;

    // Added a use foreground toggle
    useForegroundColorCheckBox = new QCheckBox(tr("Use foreground color"), this);
    connect(useForegroundColorCheckBox, SIGNAL(toggled(bool)),
            this, SLOT(useForegroundColorToggled(bool)));
    cLayout->addWidget(useForegroundColorCheckBox, row, 0, 1, 4);
    ++row;

    // Add controls for the line color.
    colorLabel = new QLabel(tr("Line color"), this);
    cLayout->addWidget(colorLabel, row, 0, Qt::AlignRight);

    colorButton = new QvisColorButton(this);
    connect(colorButton, SIGNAL(selectedColor(const QColor &)),
            this, SLOT(colorChanged(const QColor &)));
    cLayout->addWidget(colorButton, row, 1);

    // Add controls for the line opacity.
    opacitySlider = new QvisOpacitySlider(0, 255, 10, 0, this);
    connect(opacitySlider, SIGNAL(valueChanged(int)),
            this, SLOT(opacityChanged(int)));
    cLayout->addWidget(opacitySlider, row, 2, 1, 2);
    ++row;

    // Beginning arrow control.
    beginArrowComboBox = new QComboBox(this);
    beginArrowComboBox->addItem(tr("None"));
    beginArrowComboBox->addItem(tr("Line"));
    beginArrowComboBox->addItem(tr("Solid"));
    beginArrowComboBox->setEditable(false);
    connect(beginArrowComboBox, SIGNAL(activated(int)),
            this, SLOT(beginArrowChanged(int)));
    cLayout->addWidget(beginArrowComboBox, row, 1, 1, 3);
    cLayout->addWidget(new QLabel(tr("Begin arrow"), this), row, 0);
    ++row;

    // Beginning arrow control.
    endArrowComboBox = new QComboBox(this);
    endArrowComboBox->addItem(tr("None"));
    endArrowComboBox->addItem(tr("Line"));
    endArrowComboBox->addItem(tr("Solid"));
    endArrowComboBox->setEditable(false);
    connect(endArrowComboBox, SIGNAL(activated(int)),
            this, SLOT(endArrowChanged(int)));
    cLayout->addWidget(endArrowComboBox, row, 1, 1, 3);
    cLayout->addWidget(new QLabel(tr("End arrow"), this), row, 0);
    ++row;

    // Added a visibility toggle
    visibleCheckBox = new QCheckBox(tr("Visible"), this);
    connect(visibleCheckBox, SIGNAL(toggled(bool)),
            this, SLOT(visibilityToggled(bool)));
    cLayout->addWidget(visibleCheckBox, row, 0);
}

// ****************************************************************************
// Method: QvisLine2DInterface::~QvisLine2DInterface
//
// Purpose: 
//   Destructor for the QvisLine2DInterface class.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:31:23 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
QvisLine2DInterface::~QvisLine2DInterface()
{
}

// ****************************************************************************
// Method: QvisLine2DInterface::GetMenuText
//
// Purpose: 
//   Returns the text to use in the annotation list box.
//
// Arguments:
//   annot : The annotation object to use for extra information.
//
// Returns:    The text to use in the annotation list box.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:31:41 PDT 2004
//
// Modifications:
//   Brad Whitlock, Mon Jul 21 10:45:24 PDT 2008
//   Qt 4.
//
// ****************************************************************************
QString
QvisLine2DInterface::GetMenuText(const AnnotationObject &annot) const
{
    QString retval;
    if(annot.GetText().size() > 0)
        retval = QString("%1 - %2").arg(GetName()).arg(annot.GetText()[0].c_str());
    else
        retval = GetName();

    return retval;
}

// ****************************************************************************
// Method: QvisLine2DInterface::UpdateControls
//
// Purpose: 
//   Updates the controls in the interface using the data in the Annotation
//   object pointed to by the annot pointer.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:31:46 PDT 2004
//
// Modifications:
//   Brad Whitlock, Mon Jul 21 10:46:14 PDT 2008
//   Qt 4.
//
//   Kathleen Biagas, Mon Jul 13 13:09:20 PDT 2015
//   Add useForegroundcolor, colorLabel.
//
// ****************************************************************************
void
QvisLine2DInterface::UpdateControls()
{
    // Set the start position.
    positionStartEdit->setPosition(annot->GetPosition()[0],
                                   annot->GetPosition()[1]);
    
    // Set the end position.
    positionEndEdit->setPosition(annot->GetPosition2()[0],
                                 annot->GetPosition2()[1]);

    // Set the values for the width and style 
    widthWidget->blockSignals(true);
    widthWidget->SetLineWidth(annot->GetOptions().GetEntry("width")->AsInt());
    widthWidget->blockSignals(false);

    // Set the begin and end arrow styles.
    beginArrowComboBox->blockSignals(true);
    endArrowComboBox->blockSignals(true);
    beginArrowComboBox->setCurrentIndex(annot->GetOptions().GetEntry("beginArrow")->AsInt());
    endArrowComboBox->setCurrentIndex(annot->GetOptions().GetEntry("endArrow")->AsInt());
    beginArrowComboBox->blockSignals(false);
    endArrowComboBox->blockSignals(false);

    // Set the use foreground color check box.
    useForegroundColorCheckBox->blockSignals(true);
    useForegroundColorCheckBox->setChecked(annot->GetUseForegroundForTextColor());
    useForegroundColorCheckBox->blockSignals(false);

    // Change color and opacity.
    colorButton->blockSignals(true);
    opacitySlider->blockSignals(true);

    if (annot->GetUseForegroundForTextColor())
    {    
        QColor tmp(255,255,255);
        colorButton->setButtonColor(tmp);
        colorLabel->setEnabled(false);
        colorButton->setEnabled(false);
        opacitySlider->setGradientColor(tmp);
    }
    else
    {
        QColor tc(annot->GetColor1().Red(),
                  annot->GetColor1().Green(),
                  annot->GetColor1().Blue());
        colorButton->setButtonColor(tc);
        colorLabel->setEnabled(true);
        colorButton->setEnabled(true);
        opacitySlider->setGradientColor(tc);
        opacitySlider->setValue(annot->GetColor1().Alpha());
    }
    opacitySlider->blockSignals(false);
    colorButton->blockSignals(false);

    // Set the visible check box.
    visibleCheckBox->blockSignals(true);
    visibleCheckBox->setChecked(annot->GetVisible());
    visibleCheckBox->blockSignals(false);
}

// ****************************************************************************
// Method: QvisLine2DInterface::GetCurrentValues
//
// Purpose: 
//   Gets the current values for the text fields.
//
// Arguments:
//   which_widget : The widget for which we're getting the values. -1 for all.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:32:04 PDT 2004
//
// Modifications:
//   Brad Whitlock, Mon Mar 6 11:08:06 PDT 2006
//   I added code to make sure that the end points get recorded.
//
// ****************************************************************************

void
QvisLine2DInterface::GetCurrentValues(int which_widget)
{
    bool doAll = (which_widget == -1);

    if(which_widget == 0 || doAll)
    {
        // Get the new position
        GetScreenPosition(positionStartEdit, tr("Start"));
    }

    if(which_widget == 1 || doAll)
    {
        // Get the new position
        GetScreenPosition2(positionEndEdit, tr("End"));
    }

}

//
// Qt Slot functions
//

// ****************************************************************************
// Method: QvisLine2DInterface::positionChanged
//
// Purpose: 
//   This is a Qt slot function that is called when return is pressed in the 
//   position line edit.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:32:19 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
void
QvisLine2DInterface::positionStartChanged(double x, double y)
{
    double pos[] = {x, y, 0.};
    annot->SetPosition(pos);
}

// ****************************************************************************
// Method: QvisLine2DInterface::positionChanged
//
// Purpose: 
//   This is a Qt slot function that is called when return is pressed in the 
//   position line edit.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:32:25 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
void
QvisLine2DInterface::positionEndChanged(double x, double y)
{
    double pos[] = {x, y, 0.};
    annot->SetPosition2(pos);
}

// ****************************************************************************
// Method: QvisLine2DInterface::beginArrowChanged
//
// Purpose:
//   Called when the begin arrow is changed.
//
// Arguments:
//   i:    The type of arrow to use.
//
// Returns:
//
// Note:
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:32:34 PDT 2004
//
// Modifications:
//
// ****************************************************************************
void
QvisLine2DInterface::beginArrowChanged(int i)
{
    annot->GetOptions().GetEntry("beginArrow")->SetValue(i);
}

// ****************************************************************************
// Method: QvisLine2DInterface::endArrowChanged
//
// Purpose:
//   Called when the end arrow is changed.
//
// Arguments:
//   i:    The type of arrow to use.
//
// Returns:
//
// Note:
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:32:34 PDT 2004
//
// Modifications:
//
// ****************************************************************************
void
QvisLine2DInterface::endArrowChanged(int i)
{
    annot->GetOptions().GetEntry("endArrow")->SetValue(i);
}

// ****************************************************************************
// Method: QvisLine2DInterface::widthChanged
//
// Purpose: 
//   This is a Qt slot function that is called when the value of the width
//   spin box changes.
//
// Arguments:
//   w : The new width in percent.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:33:26 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
void
QvisLine2DInterface::widthChanged(int w)
{
    annot->GetOptions().GetEntry("width")->SetValue(w);
}

// ****************************************************************************
// Method: QvisLine2DInterface::colorChanged
//
// Purpose: 
//   This is a Qt slot function that is called when a new color is
//   selected.
//
// Arguments:
//   c : The new start color.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:33:35 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
void
QvisLine2DInterface::colorChanged(const QColor &c)
{
    int a = annot->GetColor1().Alpha();
    ColorAttribute tc(c.red(), c.green(), c.blue(), a);
    annot->SetColor1(tc);
}

// ****************************************************************************
// Method: QvisLine2DInterface::opacityChanged
//
// Purpose: 
//   This is a Qt slot function that is called when a new opacity is
//   selected.
//
// Arguments:
//   opacity : The new start opacity.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:33:47 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
void
QvisLine2DInterface::opacityChanged(int opacity)
{
    ColorAttribute tc(annot->GetColor1());
    tc.SetAlpha(opacity);
    annot->SetColor1(tc);
}

// ****************************************************************************
// Method: QvisLine2DInterface::visibilityToggled
//
// Purpose: 
//   This is a Qt slot function that is called when the visibility toggle is
//   changed.
//
// Arguments:
//   val : The visibility flag.
//
// Programmer: John C. Anderson
// Creation:   Fri Sep 03 09:34:03 PDT 2004
//
// Modifications:
//   
// ****************************************************************************
void
QvisLine2DInterface::visibilityToggled(bool val)
{
    annot->SetVisible(val);
}

// ****************************************************************************
// Method: QvisLine2DInterface::useForegroundColorToggled
//
// Purpose: 
//   This is a Qt slot function that is called when the useForegroundColor
//   check box is clicked.
//
// Arguments:
//   val : The new setting for useForegroundColor
//
// Programmer: Kathleen Biagas 
// Creation:   July 13, 2015
//
// Modifications:
//   
// ****************************************************************************

void
QvisLine2DInterface::useForegroundColorToggled(bool val)
{
    annot->SetUseForegroundForTextColor(val);
    Apply();
}

