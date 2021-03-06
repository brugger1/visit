// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef QVISSEEDMEWINDOW_H
#define QVISSEEDMEWINDOW_H

#include <SeedMeAttributes.h>
#include <QvisPostableWindowSimpleObserver.h>

class QButtonGroup;
class QCheckBox;
class QLabel;
class QTextBrowser;
class QLineEdit;
class QTabWidget;
class QPushButton;
class QFileSystemWatcher;
class QComboBox;

// ****************************************************************************
// Class: QvisSeedMeWindow
//
// Purpose:
//    Defines QvisSeedMeWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//    Kathleen Biagas, Mon Aug 14 10:43:08 MST 2017
//    Added clearTabsOnClose, hide slot, closeEvent.
//
// ****************************************************************************

class QvisSeedMeWindow : public QvisPostableWindowSimpleObserver
{
    Q_OBJECT
public:
    QvisSeedMeWindow(SeedMeAttributes *subj,
                     const QString &caption = QString::null,
                     const QString &shortName = QString::null,
                     QvisNotepadArea *notepad = 0);
    virtual ~QvisSeedMeWindow();
    virtual void CreateWindowContents();

    virtual void SubjectRemoved(Subject *subj);
signals:
    void runCommand(const QString &);
public slots:
    virtual void apply();
    virtual void hide();
protected:
    void closeEvent(QCloseEvent *event);
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);
private slots:
    void collectionModeChanged(int val);
    void collectionIDProcessText();
    void sharingChanged(int val);
    void collectionTitleProcessText();
    void collectionDescriptionProcessText();
    void overwriteFilesChanged(bool val);
    void keyValueProcessText();
    void collectionEmailsProcessText();
    void uploadSequenceFileChanged(bool val);
    void sequenceTitleProcessText();
    void sequenceDescriptionProcessText();
    void createVideoChanged(bool val);
    void frameRateProcessText();
    void browse();
    void browseApiKey();
    void queryActionChanged(int val);
    void queryColIDProcessText();
    void queryKeyValueProcessText();
    void queryCollectionValuesChanged(int val);
    void downloadCollectionIDProcessText();
    void downloadTypeChanged(int val);
    void downloadNameProcessText();
    void directoryChanged(const QString & path);
    void quickSharingChanged(int val);
    void quickCollectionTitleProcessText();
    void quickCollectionEmailsProcessText();
    void quickFrameRateProcessText();
    void quickDownloadTypeChanged(int val);
    void quickDownload();
    void clearTabsOnCloseChanged(bool val);


    void ResetForm();
    void ClearLog();
private:
    void ResetForm2(int);
    QWidget *CreateQuickUploadTab();
    QWidget *CreateUploadTab();
    QWidget *CreateQueryTab();
    QWidget *CreateDownloadTab();
    QWidget *CreateSettingsTab();

    QTabWidget   *tabs;

    // Upload tab
    QWidget      *collectionMode;
    QButtonGroup *collectionModeButtonGroup;
    QLineEdit *collectionID;
    QWidget      *sharing;
    QButtonGroup *sharingButtonGroup;
    QLineEdit *collectionTitle;
    QLineEdit *collectionDescription;
    QLineEdit *collectionCredits;
    QLineEdit *collectionLicense;
    QLineEdit *keyValue;
    QLineEdit *collectionEmails;
    QCheckBox *overwriteFiles;
    QCheckBox *uploadSequenceFile;
    QLineEdit *sequenceTitle;
    QLineEdit *sequenceDescription;
    QCheckBox *createVideo;
    QLineEdit *frameRate;
    QWidget      *queryAction;
    QButtonGroup *queryActionButtonGroup;
    QLineEdit *queryColID;
    QLineEdit *queryKeyValue;
    QWidget      *queryCollectionValues;
    QButtonGroup *queryCollectionValuesButtonGroup;
    QLineEdit *downloadCollectionID;
    QWidget      *downloadType;
    QButtonGroup *downloadTypeButtonGroup;
    QLineEdit *downloadName;
    QWidget      *quickSharing;
    QButtonGroup *quickSharingButtonGroup;
    QLineEdit *quickCollectionTitle;
    QLineEdit *quickCollectionEmails;
    QLineEdit *quickFrameRate;
    QComboBox      *quickDownloadType;
    QLabel *collectionModeLabel;
    QLabel *collectionIDLabel;
    QLabel *sharingLabel;
    QLabel *collectionTitleLabel;
    QLabel *collectionDescriptionLabel;
    QLabel *keyValueLabel;
    QLabel *collectionCreditsLabel;
    QLabel *collectionLicenseLabel;
    QLabel *collectionEmailsLabel;
    QLabel *currentTitleLabel;
    QLabel *currentDescriptionLabel;
    QLabel *sequenceTitleLabel;
    QLabel *sequenceDescriptionLabel;
    QLabel *frameRateLabel;
    QLabel *queryActionLabel;
    QLabel *queryColIDLabel;
    QLabel *queryKeyValueLabel;
    QLabel *queryCollectionValuesLabel;
    QLabel *downloadCollectionIDLabel;
    QLabel *downloadTypeLabel;
    QLabel *downloadNameLabel;

    QLabel *collectionsLink, *quickCollectionsLink;

    QLabel *helpLabel;
    QLabel *helpLabelWarning;

    QLabel *confLocationLabel;

    QFileSystemWatcher *seedmeWatcher;

    QTextBrowser *statusLabel;

    SeedMeAttributes *atts;

    QStringList uploadFiles;

    QPushButton *resetFormButton, *clearLogButton, *quickDownloadButton, *submitButton;

    QCheckBox *clearTabsOnClose;

    void updateStatus(QString str);

    QString apikeyFile;

    const char* configFile;
};


#endif
