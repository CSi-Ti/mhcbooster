# pip install PySide2
import os
import re
import subprocess
import sys
import webbrowser
import logging
import time
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QPixmap, QIcon, QTextCursor
from PySide6.QtWidgets import (QApplication, QWidget, QPushButton, QLabel,
                               QVBoxLayout, QLineEdit, QFileDialog, QHBoxLayout,
                               QCheckBox, QGridLayout, QSpinBox, QGroupBox,
                               QMessageBox, QTextEdit, QTabWidget, QStackedWidget, QSizePolicy, QRadioButton,
                               QDoubleSpinBox, QProgressBar)

from src.utils.package_installer import *

ROOT_PATH = Path(__file__).parent.parent.parent
sys.path.append(ROOT_PATH.as_posix())
from src import __version__


def grid_layout(label, elements, n_same_row=4):
    g_layout = QGridLayout()
    g_layout.setHorizontalSpacing(10)
    g_layout.setVerticalSpacing(3)
    for i, checkbox in enumerate(elements):
        row = i // n_same_row  # Every 5 checkboxes will be placed in a new row
        col = i % n_same_row  # Columns will repeat after every 5 checkboxes (like a 5-column grid)
        g_layout.addWidget(checkbox, row, col)
        g_layout.setColumnMinimumWidth(col, 220)
    h_layout = QHBoxLayout()
    h_layout.addWidget(label)
    h_layout.setAlignment(label, Qt.AlignTop)
    h_layout.addLayout(g_layout)
    h_layout.setAlignment(Qt.AlignLeft)
    return h_layout


class MhcBoosterGUI(QWidget):
    def __init__(self):
        super().__init__()

        self.setStyleSheet("""
                QGroupBox {
                    border: 1px solid #A9A9A9;
                    margin-top: 1ex;
                    padding: 5px;
                    font: bold 12px;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    subcontrol-position: top left; /* position at the top left edge */
                    padding: 0 3px; /* padding from the border */
                    left: 10px;
                }
                QTabWidget::pane {
                    border: 1px solid lightgray;  /* Remove the tab box */
                    border-left: none;
                    border-right: none;
                    border-bottom: none;
                }
            """)

        # GUI window
        self.setWindowTitle('MhcBooster')
        self.setWindowIcon(QIcon(str(Path(__file__).parent/'caronlab_icon.png')))
        # self.setGeometry(100, 100, 800, 600)
        layout = QVBoxLayout()
        # layout.setContentsMargins(50, 20, 50, 10) # left, top, right, bottom
        layout.setSpacing(20)

        ### INTRODUCTION
        logo_lab_label = QLabel()
        logo_pix_map = QPixmap(str(Path(__file__).parent/'caronlab.png')).scaled(200, 150, Qt.KeepAspectRatio)
        logo_lab_label.setPixmap(logo_pix_map)
        logo_lab_label.resize(logo_pix_map.size())
        intro_label = QLabel('The Introduction of MhcBooster should be here. GitHub. Tutorial. Cite. CaronLab.')
        intro_layout = QHBoxLayout()
        intro_layout.addWidget(logo_lab_label)
        intro_layout.addWidget(intro_label)
        layout.addLayout(intro_layout)


        self.tab_widget = QTabWidget()
        self.add_main_tab()
        self.add_reporter_tab()
        self.add_config_tab()
        layout.addWidget(self.tab_widget)

        ### Footnote
        foot_label = QLabel('CaronLab 2024')
        foot_label.setAlignment(Qt.AlignRight)
        layout.addWidget(foot_label)

        self.setLayout(layout)


    def add_main_tab(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(30, 20, 30, 20) # left, top, right, bottom
        main_layout.setSpacing(20)

        ### FILE MANAGEMENT
        file_groupbox = QGroupBox('Input / Output')
        file_group_layout = QVBoxLayout()

        input_format_layout = QHBoxLayout()
        input_format_layout.setContentsMargins(0, 0, 0, 0)
        input_format_layout.setSpacing(30)
        self.pin_radiobutton = QRadioButton('Run from PSM')
        self.raw_radiobutton = QRadioButton('Run from RAW (MSFragger)')
        self.pin_radiobutton.setChecked(True)
        input_format_layout.addWidget(self.pin_radiobutton)
        input_format_layout.addWidget(self.raw_radiobutton)
        input_format_layout.setAlignment(Qt.AlignLeft)
        file_group_layout.addLayout(input_format_layout)

        psm_group_layout = QVBoxLayout()
        psm_group_layout.setContentsMargins(0, 0, 0, 0)
        psm_group_layout.setSpacing(5)

        # PIN folder
        psm_path_label = QLabel('PSM folder: \t')
        self.psm_inputbox = QLineEdit()
        self.psm_inputbox.setPlaceholderText("Select input folder containing .pin files from Comet or MSFragger...")
        self.psm_button = QPushButton("Select")
        self.psm_button.clicked.connect(self.open_folder_dialog)
        psm_layout = QHBoxLayout()
        psm_layout.addWidget(psm_path_label)
        psm_layout.addWidget(self.psm_inputbox)
        psm_layout.addWidget(self.psm_button)
        psm_group_layout.addLayout(psm_layout)

        # MzML folder
        mzml_label = QLabel('mzML folder: \t')
        self.mzml_inputbox = QLineEdit()
        self.mzml_inputbox.setPlaceholderText('Select mzML folder containing .mzML files with the same name as PSM files...')
        self.mzml_button = QPushButton("Select")
        self.mzml_button.clicked.connect(self.open_folder_dialog)
        mzml_layout = QHBoxLayout()
        mzml_layout.addWidget(mzml_label)
        mzml_layout.addWidget(self.mzml_inputbox)
        mzml_layout.addWidget(self.mzml_button)
        psm_group_layout.addLayout(mzml_layout)

        # Output folder
        psm_output_label = QLabel('Output folder: \t')
        self.psm_output_inputbox = QLineEdit()
        self.psm_output_inputbox.setPlaceholderText('Select output folder...')
        self.psm_output_button = QPushButton("Select")
        self.psm_output_button.clicked.connect(self.open_folder_dialog)
        psm_output_layout = QHBoxLayout()
        psm_output_layout.addWidget(psm_output_label)
        psm_output_layout.addWidget(self.psm_output_inputbox)
        psm_output_layout.addWidget(self.psm_output_button)
        psm_group_layout.addLayout(psm_output_layout)

        raw_group_layout = QVBoxLayout()
        raw_group_layout.setContentsMargins(0, 0, 0, 0)
        raw_group_layout.setSpacing(5)

        # RAW folder
        raw_path_label = QLabel('RAW folder: \t')
        self.raw_inputbox = QLineEdit()
        self.raw_inputbox.setPlaceholderText("Select input folder containing .raw/.d files...")
        self.raw_button = QPushButton("Select")
        self.raw_button.clicked.connect(self.open_folder_dialog)
        raw_layout = QHBoxLayout()
        raw_layout.addWidget(raw_path_label)
        raw_layout.addWidget(self.raw_inputbox)
        raw_layout.addWidget(self.raw_button)
        raw_group_layout.addLayout(raw_layout)

        # fasta file
        fasta_label = QLabel('FASTA file: \t')
        self.fasta_inputbox = QLineEdit()
        self.fasta_inputbox.setPlaceholderText('Select .fasta or .fasta.fas file...')
        self.fasta_button = QPushButton("Select")
        self.fasta_button.clicked.connect(self.open_file_dialog)
        fasta_layout = QHBoxLayout()
        fasta_layout.addWidget(fasta_label)
        fasta_layout.addWidget(self.fasta_inputbox)
        fasta_layout.addWidget(self.fasta_button)
        raw_group_layout.addLayout(fasta_layout)

        # parameter file
        param_label = QLabel('Parameter file: \t')
        self.param_inputbox = QLineEdit()
        self.param_inputbox.setPlaceholderText('Select fragger.params file...')
        self.param_button = QPushButton("Select")
        self.param_button.clicked.connect(self.open_file_dialog)
        param_layout = QHBoxLayout()
        param_layout.addWidget(param_label)
        param_layout.addWidget(self.param_inputbox)
        param_layout.addWidget(self.param_button)
        raw_group_layout.addLayout(param_layout)

        # Output folder
        raw_output_label = QLabel('Output folder: \t')
        self.raw_output_inputbox = QLineEdit()
        self.raw_output_inputbox.setPlaceholderText('Select output folder...')
        self.raw_output_button = QPushButton("Select")
        self.raw_output_button.clicked.connect(self.open_folder_dialog)
        raw_output_layout = QHBoxLayout()
        raw_output_layout.addWidget(raw_output_label)
        raw_output_layout.addWidget(self.raw_output_inputbox)
        raw_output_layout.addWidget(self.raw_output_button)
        raw_group_layout.addLayout(raw_output_layout)

        pin_group_widget = QWidget()
        pin_group_widget.setLayout(psm_group_layout)
        raw_group_widget = QWidget()
        raw_group_widget.setLayout(raw_group_layout)

        self.input_stacked_widget = QStackedWidget()
        self.input_stacked_widget.addWidget(pin_group_widget)
        self.input_stacked_widget.addWidget(raw_group_widget)

        file_group_layout.addWidget(self.input_stacked_widget)
        self.pin_radiobutton.clicked.connect(lambda: self.input_stacked_widget.setCurrentIndex(0))
        self.raw_radiobutton.clicked.connect(lambda: self.input_stacked_widget.setCurrentIndex(1))
        file_groupbox.setLayout(file_group_layout)
        main_layout.addWidget(file_groupbox)

        ### MHC specific SCORES
        mhc_groupbox = QGroupBox('MHC Predictors')
        mhc_group_layout = QVBoxLayout()
        mhc_group_layout.insertSpacing(0, 5)
        mhc_group_layout.setSpacing(10)

        # APP score
        mhc_I_label = QLabel('MHC-I Score:\t')
        mhc_I_models = ['NetMHCpan', 'MHCflurry', 'BigMHC']
        self.checkboxes_mhc_I = [QCheckBox(model) for model in mhc_I_models]
        for checkbox in self.checkboxes_mhc_I:
            checkbox.toggled.connect(self.on_mhc_I_checkbox_toggled)
        mhc_I_layout = grid_layout(mhc_I_label, self.checkboxes_mhc_I)
        mhc_group_layout.addLayout(mhc_I_layout)
        mhc_II_label = QLabel('MHC-II Score:\t')
        mhc_II_models = ['NetMHCIIpan', 'MixMHC2pred']
        self.checkboxes_mhc_II = [QCheckBox(model) for model in mhc_II_models]
        for checkbox in self.checkboxes_mhc_II:
            checkbox.toggled.connect(self.on_mhc_II_checkbox_toggled)
        mhc_II_layout = grid_layout(mhc_II_label, self.checkboxes_mhc_II)
        mhc_group_layout.addLayout(mhc_II_layout)

        # Alleles
        allele_label = QLabel('Alleles: \t   ')
        self.allele_inputbox = QLineEdit()
        self.allele_inputbox.setPlaceholderText('Input alleles (e.g. HLA-A0101; DQB1*05:01) or Select allele map file...')
        self.allele_button = QPushButton("Select")
        self.allele_button.clicked.connect(self.open_file_dialog)
        allele_layout = QHBoxLayout()
        allele_layout.addWidget(allele_label)
        allele_layout.addWidget(self.allele_inputbox)
        allele_layout.addWidget(self.allele_button)
        mhc_group_layout.addLayout(allele_layout)

        mhc_groupbox.setLayout(mhc_group_layout)
        main_layout.addWidget(mhc_groupbox)

        ### GENERAL SCORES
        gs_groupbox = QGroupBox('General Predictors')
        gs_group_layout = QVBoxLayout()
        gs_group_layout.insertSpacing(0, 5)
        gs_group_layout.setSpacing(20)

        # RT score
        rt_label = QLabel('RT Score: \t')
        rt_models = ['AutoRT', 'Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Chronologer_RT',
                     'Prosit_2019_irt', 'Prosit_2024_irt_cit', 'Prosit_2020_irt_TMT']
        self.checkboxes_rt = [QCheckBox(model) for model in rt_models]
        rt_layout = grid_layout(rt_label, self.checkboxes_rt)
        gs_group_layout.addLayout(rt_layout)

        # MS2 score
        ms2_label = QLabel('MS2 Score:\t')
        unsuitable_ms2_models = ['UniSpec', 'Prosit_2024_intensity_XL_NMS2', 'Prosit_2023_intensity_XL_CMS2',
                                 'Prosit_2023_intensity_XL_CMS3']
        ms2_models = ['AlphaPeptDeep_ms2_generic', 'ms2pip_HCD2021', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2023',
                      'ms2pip_timsTOF2024', 'ms2pip_iTRAQphospho', 'ms2pip_TTOF5600', 'ms2pip_CID_TMT',
                      'Prosit_2019_intensity', 'Prosit_2020_intensity_HCD', 'Prosit_2020_intensity_CID',
                      'Prosit_2023_intensity_timsTOF', 'Prosit_2024_intensity_cit', 'Prosit_2020_intensity_TMT']

        self.checkboxes_ms2 = [QCheckBox(model) for model in ms2_models]
        self.checkboxes_ms2.insert(7, QLabel(''))
        ms2_layout = grid_layout(ms2_label, self.checkboxes_ms2)
        gs_group_layout.addLayout(ms2_layout)

        # CCS score
        ccs_label = QLabel('CCS Score:\t')
        ccs_models = ['IM2Deep', 'AlphaPeptDeep_ccs_generic']
        self.checkboxes_ccs = [QCheckBox(model) for model in ccs_models]
        ccs_layout = grid_layout(ccs_label, self.checkboxes_ccs)
        gs_group_layout.addLayout(ccs_layout)

        # Peptide encoding
        pe_label = QLabel('Sequence Feature:\t')
        pe_models = ['Peptide Encoding']
        self.checkboxes_pe = [QCheckBox(model) for model in pe_models]
        pe_layout = grid_layout(pe_label, self.checkboxes_pe)
        gs_group_layout.addLayout(pe_layout)

        # Auto select
        ap_layout = QHBoxLayout()
        self.ap_checkbox = QCheckBox('Auto-predict best combination')
        self.ap_checkbox.toggled.connect(self.on_autopred_checkbox_toggled)
        ap_layout.addWidget(self.ap_checkbox)
        gs_group_layout.addLayout(ap_layout)

        gs_groupbox.setLayout(gs_group_layout)
        main_layout.addWidget(gs_groupbox)


        ### RUN PARAMS
        rp_groupbox = QGroupBox('Run Parameters')
        rp_group_layout = QVBoxLayout()
        rp_group_layout.insertSpacing(0, 5)

        p1_layout = QHBoxLayout()
        p1_layout.setAlignment(Qt.AlignLeft)

        # Fine tune
        self.ft_checkbox = QCheckBox('Fine tune')
        p1_layout.addWidget(self.ft_checkbox)
        p1_layout.setAlignment(Qt.AlignLeft)
        p1_layout.addSpacing(30)

        # Peptide length
        self.pl_checkbox = QCheckBox('Filter by length')
        self.pl_min = QSpinBox()
        self.pl_min.setRange(7, 20)
        self.pl_min.setValue(8)
        self.pl_max = QSpinBox()
        self.pl_max.setRange(7, 20)
        self.pl_max.setValue(15)
        p1_layout.addWidget(self.pl_checkbox)
        p1_layout.addWidget(self.pl_min)
        p1_layout.addWidget(QLabel('-'))
        p1_layout.addWidget(self.pl_max)
        p1_layout.addSpacing(30)

        # Koina
        koina_label = QLabel('Koina server URL: ')
        self.koina_inputbox = QLineEdit('koina.wilhelmlab.org:443')
        p1_layout.addWidget(koina_label)
        p1_layout.addWidget(self.koina_inputbox)
        p1_layout.addSpacing(30)

        # Max thread
        self.thread_label = QLabel('Threads: ')
        self.spinbox_thread = QSpinBox()
        self.spinbox_thread.setRange(1, os.cpu_count() - 1)
        self.spinbox_thread.setValue(os.cpu_count() - 1)
        p1_layout.addWidget(self.thread_label)
        p1_layout.addWidget(self.spinbox_thread)
        p1_layout.addSpacing(80)
        rp_group_layout.addLayout(p1_layout)

        rp_groupbox.setLayout(rp_group_layout)
        main_layout.addWidget(rp_groupbox)

        ### Logger
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setFixedHeight(150)
        main_layout.addWidget(self.log_output)

        ### Execution
        self.button_run = QPushButton("RUN")
        self.button_run.clicked.connect(self.on_exec_clicked)
        main_layout.addWidget(self.button_run)

        self.worker_thread = MhcBoosterWorker(commands=None)
        self.worker_thread.message.connect(self.add_log)
        self.worker_thread.finished.connect(self.worker_stop)

        main_tab = QWidget()
        main_tab.setLayout(main_layout)
        self.tab_widget.addTab(main_tab, 'Identification')


    def add_config_tab(self):
        config_layout = QVBoxLayout()
        config_layout.setContentsMargins(30, 20, 30, 20) # left, top, right, bottom
        config_layout.setSpacing(20)
        config_layout.setAlignment(Qt.AlignTop)

        third_party_groupbox = QGroupBox('Third-party tools')
        third_party_layout = QVBoxLayout()

        # Introduction
        introduction_label = QLabel('MHCBooster utilizes a variety of tools for RT, MS2, and CCS scoring.'
                                    ' Some of these tools are governed by strict licenses and must be manually'
                                    ' downloaded and installed. Please input the paths to the downloaded'
                                    ' zip files. And they will be automatically installed by pressing'
                                    ' the \'Install to MHCBooster\' button.')
        introduction_label.setWordWrap(True)
        introduction_label_layout = QHBoxLayout()
        introduction_label_layout.setContentsMargins(0, 10, 0, 20)
        introduction_label_layout.addWidget(introduction_label)
        third_party_layout.addLayout(introduction_label_layout)

        # MSFragger
        msfragger_label = QLabel('MSFragger path: \t')
        self.msfragger_inputbox = QLineEdit()
        self.msfragger_inputbox.setPlaceholderText("Select the path to MSFragger-4.1.zip ...")
        self.msfragger_browse_button = QPushButton("Browse")
        self.msfragger_browse_button.clicked.connect(self.open_file_dialog)
        self.msfragger_download_button = QPushButton("Download")
        self.msfragger_download_button.clicked.connect(lambda: webbrowser.open('https://msfragger-upgrader.nesvilab.org/upgrader/'))
        msfragger_layout = QHBoxLayout()
        msfragger_layout.addWidget(msfragger_label)
        msfragger_layout.addWidget(self.msfragger_inputbox)
        msfragger_layout.addWidget(self.msfragger_browse_button)
        msfragger_layout.addWidget(self.msfragger_download_button)
        third_party_layout.addLayout(msfragger_layout)

        # AutoRT
        autort_label = QLabel('AutoRT path: \t')
        self.autort_inputbox = QLineEdit()
        self.autort_inputbox.setPlaceholderText("Select the path to AutoRT-master.zip ...")
        self.autort_browse_button = QPushButton("Browse")
        self.autort_browse_button.clicked.connect(self.open_file_dialog)
        self.autort_download_button = QPushButton("Download")
        self.autort_download_button.clicked.connect(lambda: webbrowser.open('https://github.com/bzhanglab/AutoRT/archive/refs/heads/master.zip'))
        autort_layout = QHBoxLayout()
        autort_layout.addWidget(autort_label)
        autort_layout.addWidget(self.autort_inputbox)
        autort_layout.addWidget(self.autort_browse_button)
        autort_layout.addWidget(self.autort_download_button)
        third_party_layout.addLayout(autort_layout)

        # BigMHC
        bigmhc_label = QLabel('BigMHC path: \t')
        self.bigmhc_inputbox = QLineEdit()
        self.bigmhc_inputbox.setPlaceholderText("Select the path to bigmhc-master.zip ...")
        self.bigmhc_browse_button = QPushButton("Browse")
        self.bigmhc_browse_button.clicked.connect(self.open_file_dialog)
        self.bigmhc_download_button = QPushButton("Download")
        self.bigmhc_download_button.clicked.connect(lambda: webbrowser.open('https://github.com/KarchinLab/bigmhc/archive/refs/heads/master.zip'))
        bigmhc_layout = QHBoxLayout()
        bigmhc_layout.addWidget(bigmhc_label)
        bigmhc_layout.addWidget(self.bigmhc_inputbox)
        bigmhc_layout.addWidget(self.bigmhc_browse_button)
        bigmhc_layout.addWidget(self.bigmhc_download_button)
        third_party_layout.addLayout(bigmhc_layout)

        # NetMHCpan
        netmhcpan_label = QLabel('NetMHCpan path: \t')
        self.netmhcpan_inputbox = QLineEdit()
        self.netmhcpan_inputbox.setPlaceholderText("Select the path to netMHCpan-4.1b.Linux.tar.gz ...")
        self.netmhcpan_browse_button = QPushButton("Browse")
        self.netmhcpan_browse_button.clicked.connect(self.open_file_dialog)
        self.netmhcpan_download_button = QPushButton("Download")
        self.netmhcpan_download_button.clicked.connect(lambda: webbrowser.open('https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCpan&version=4.1&packageversion=4.1b&platform=Linux'))
        netmhcpan_layout = QHBoxLayout()
        netmhcpan_layout.addWidget(netmhcpan_label)
        netmhcpan_layout.addWidget(self.netmhcpan_inputbox)
        netmhcpan_layout.addWidget(self.netmhcpan_browse_button)
        netmhcpan_layout.addWidget(self.netmhcpan_download_button)
        third_party_layout.addLayout(netmhcpan_layout)

        # NetMHCIIpan
        netmhcIIpan_label = QLabel('NetMHCIIpan path:\t')
        self.netmhcIIpan_inputbox = QLineEdit()
        self.netmhcIIpan_inputbox.setPlaceholderText("Select the path to netMHCIIpan-4.3e.Linux.tar.gz ...")
        self.netmhcIIpan_browse_button = QPushButton("Browse")
        self.netmhcIIpan_browse_button.clicked.connect(self.open_file_dialog)
        self.netmhcIIpan_download_button = QPushButton("Download")
        self.netmhcIIpan_download_button.clicked.connect(lambda: webbrowser.open('https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCIIpan&version=4.3&packageversion=4.3e&platform=Linux'))
        netmhcIIpan_layout = QHBoxLayout()
        netmhcIIpan_layout.addWidget(netmhcIIpan_label)
        netmhcIIpan_layout.addWidget(self.netmhcIIpan_inputbox)
        netmhcIIpan_layout.addWidget(self.netmhcIIpan_browse_button)
        netmhcIIpan_layout.addWidget(self.netmhcIIpan_download_button)
        third_party_layout.addLayout(netmhcIIpan_layout)

        # MixMHC2pred
        mixmhc2pred_label = QLabel('MixMHC2pred path:')
        self.mixmhc2pred_inputbox = QLineEdit()
        self.mixmhc2pred_inputbox.setPlaceholderText("Select the path to MixMHC2pred-2.0.zip ...")
        self.mixmhc2pred_browse_button = QPushButton("Browse")
        self.mixmhc2pred_browse_button.clicked.connect(self.open_file_dialog)
        self.mixmhc2pred_download_button = QPushButton("Download")
        self.mixmhc2pred_download_button.clicked.connect(lambda: webbrowser.open('https://github.com/GfellerLab/MixMHC2pred/releases/download/v2.0.2.2/MixMHC2pred-2.0.zip'))
        mixmhc2pred_layout = QHBoxLayout()
        mixmhc2pred_layout.addWidget(mixmhc2pred_label)
        mixmhc2pred_layout.addWidget(self.mixmhc2pred_inputbox)
        mixmhc2pred_layout.addWidget(self.mixmhc2pred_browse_button)
        mixmhc2pred_layout.addWidget(self.mixmhc2pred_download_button)
        third_party_layout.addLayout(mixmhc2pred_layout)
        
        extract_layout = QHBoxLayout()
        extract_layout.setContentsMargins(0, 10, 0, 0)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 6)  # Set range to 0 to make it indeterminate
        self.progress_bar.setTextVisible(False)  # Hide the text inside the progress bar
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(False)
        self.extract_button = QPushButton("Install to MHCBooster")
        self.extract_button.clicked.connect(self.on_install_clicked)
        self.refresh_third_party_status()
        extract_layout.addWidget(self.progress_bar)
        extract_layout.addWidget(self.extract_button)
        extract_layout.setAlignment(Qt.AlignRight)
        third_party_layout.addLayout(extract_layout)

        third_party_groupbox.setLayout(third_party_layout)
        config_layout.addWidget(third_party_groupbox)

        license_groupbox = QGroupBox('License')
        license_layout = QVBoxLayout()
        license_text = QLabel('MHCBooster is an open-source software tool released under the GNU General Public'
                              ' License (GPL) version 3. This means that you are free to use, modify, and distribute'
                              ' the software, as long as you adhere to the terms and conditions set forth by the GPL-3'
                              ' license. Additionally, when using MHCBooster, please ensure that you comply with the'
                              ' licenses of any third-party tools or libraries integrated with the software. These'
                              ' third-party components may be subject to different licensing agreements, and it is'
                              ' your responsibility to review and follow the relevant terms for each of them. By using'
                              ' MHCBooster, you agree to abide by the obligations of both the GPL-3 and any applicable'
                              ' third-party licenses.')
        license_text.setWordWrap(True)
        license_layout.addWidget(license_text)
        license_groupbox.setLayout(license_layout)
        config_layout.addWidget(license_groupbox)

        cite_groupbox = QGroupBox('How to cite')
        cite_layout = QVBoxLayout()
        cite = ('MHCBooster: <br>'
                'cite info<br>'
                'Third-party tools:<br>'
                'MSFragger:<br>'
                'Koina:<br>')
        cite_text = QTextEdit()
        cite_text.setHtml(cite)
        cite_text.setReadOnly(True)
        cite_layout.addWidget(cite_text)
        cite_groupbox.setLayout(cite_layout)
        config_layout.addWidget(cite_groupbox)

        config_tab = QWidget()
        config_tab.setLayout(config_layout)
        self.tab_widget.addTab(config_tab, 'Configuration')


    def add_reporter_tab(self):

        reporter_layout = QVBoxLayout()
        reporter_layout.setContentsMargins(30, 20, 30, 20) # left, top, right, bottom
        reporter_layout.setSpacing(20)
        reporter_layout.setAlignment(Qt.AlignTop)

        # # Introduction
        # introduction_label = QLabel('Add some descriptions and links here')
        # introduction_label.setWordWrap(True)
        # introduction_label_layout = QHBoxLayout()
        # introduction_label_layout.setContentsMargins(0, 10, 0, 20)
        # introduction_label_layout.addWidget(introduction_label)
        # reporter_layout.addLayout(introduction_label_layout)


        ### Input / Output
        input_output_groupbox = QGroupBox('Input / Output')
        input_output_layout = QVBoxLayout()
        input_output_layout.insertSpacing(0, 5)

        # Input result folder (support MHCBooster, Percolator, Mokapot, Philosopher?)
        result_path_label = QLabel('Result folder: \t')
        self.result_inputbox = QLineEdit()
        self.result_inputbox.setPlaceholderText("Select the MHCBooster results folder for which you want a custom report...")
        self.result_button = QPushButton("Select")
        self.result_button.clicked.connect(self.open_folder_dialog)
        result_layout = QHBoxLayout()
        result_layout.addWidget(result_path_label)
        result_layout.addWidget(self.result_inputbox)
        result_layout.addWidget(self.result_button)
        input_output_layout.addLayout(result_layout)
        input_output_groupbox.setLayout(input_output_layout)
        reporter_layout.addWidget(input_output_groupbox)

        ### Peptide
        # Alignment & MBR?
        # Missing value imputation？
        # Intensity normalization?
        pept_info_groupbox = QGroupBox('Peptide Info')
        pept_info_layout = QVBoxLayout()
        pept_info_layout.insertSpacing(0, 5)
        pept_info_layout.setSpacing(20)

        self.binder_checkbox = QCheckBox('Binder prediction')

        self.netmhcpan_checkbox = QCheckBox('NetMHCpan')
        netmhcpan_strong_label = QLabel('\tStrong Rank (%)')
        self.netmhcpan_strong_spinbox = QDoubleSpinBox()
        self.netmhcpan_strong_spinbox.setRange(0.001, 100)
        self.netmhcpan_strong_spinbox.setValue(0.5)
        self.netmhcpan_strong_spinbox.setSingleStep(0.1)
        netmhcpan_weak_label = QLabel('Weak Rank (%)')
        self.netmhcpan_weak_spinbox = QDoubleSpinBox()
        self.netmhcpan_weak_spinbox.setRange(0.001, 100)
        self.netmhcpan_weak_spinbox.setValue(2)
        self.netmhcpan_weak_spinbox.setSingleStep(0.5)

        self.mhcflurry_checkbox = QCheckBox('MHCflurry')
        mhcflurry_strong_label = QLabel('\tStrong Rank (%)')
        self.mhcflurry_strong_spinbox = QDoubleSpinBox()
        self.mhcflurry_strong_spinbox.setRange(0.001, 100)
        self.mhcflurry_strong_spinbox.setValue(0.5)
        self.mhcflurry_strong_spinbox.setSingleStep(0.1)
        mhcflurry_weak_label = QLabel('Weak Rank (%)')
        self.mhcflurry_weak_spinbox = QDoubleSpinBox()
        self.mhcflurry_weak_spinbox.setRange(0.001, 100)
        self.mhcflurry_weak_spinbox.setValue(2)
        self.mhcflurry_weak_spinbox.setSingleStep(0.5)

        self.bigmhc_checkbox = QCheckBox('BigMHC')
        bigmhc_strong_label = QLabel('\tStrong Rank (%)')
        self.bigmhc_strong_spinbox = QDoubleSpinBox()
        self.bigmhc_strong_spinbox.setRange(0.001, 100)
        self.bigmhc_strong_spinbox.setValue(0.5)
        self.bigmhc_strong_spinbox.setSingleStep(0.1)
        bigmhc_weak_label = QLabel('Weak Rank (%)')
        self.bigmhc_weak_spinbox = QDoubleSpinBox()
        self.bigmhc_weak_spinbox.setRange(0.001, 100)
        self.bigmhc_weak_spinbox.setValue(2)
        self.bigmhc_weak_spinbox.setSingleStep(0.5)

        self.netmhcIIpan_checkbox = QCheckBox('NetMHCIIpan')
        netmhcIIpan_strong_label = QLabel('\tStrong Rank (%)')
        self.netmhcIIpan_strong_spinbox = QDoubleSpinBox()
        self.netmhcIIpan_strong_spinbox.setRange(0.001, 100)
        self.netmhcIIpan_strong_spinbox.setValue(0.5)
        self.netmhcIIpan_strong_spinbox.setSingleStep(0.1)
        netmhcIIpan_weak_label = QLabel('Weak Rank (%)')
        self.netmhcIIpan_weak_spinbox = QDoubleSpinBox()
        self.netmhcIIpan_weak_spinbox.setRange(0.001, 100)
        self.netmhcIIpan_weak_spinbox.setValue(2)
        self.netmhcIIpan_weak_spinbox.setSingleStep(0.5)

        self.mixmhc2pred_checkbox = QCheckBox('MixMHC2pred')
        mixmhc2pred_strong_label = QLabel('\tStrong Rank (%)')
        self.mixmhc2pred_strong_spinbox = QDoubleSpinBox()
        self.mixmhc2pred_strong_spinbox.setRange(0.001, 100)
        self.mixmhc2pred_strong_spinbox.setValue(0.5)
        self.mixmhc2pred_strong_spinbox.setSingleStep(0.1)
        mixmhc2pred_weak_label = QLabel('Weak Rank (%)')
        self.mixmhc2pred_weak_spinbox = QDoubleSpinBox()
        self.mixmhc2pred_weak_spinbox.setRange(0.001, 100)
        self.mixmhc2pred_weak_spinbox.setValue(2)
        self.mixmhc2pred_weak_spinbox.setSingleStep(0.5)

        binder_layout = QHBoxLayout()
        binder_layout.setAlignment(Qt.AlignLeft)
        binder_checkbox_layout = QVBoxLayout()
        binder_checkbox_layout.setAlignment(Qt.AlignTop)
        binder_checkbox_layout.addWidget(self.binder_checkbox)
        binder_layout.addLayout(binder_checkbox_layout)
        binder_layout.addSpacing(100)

        binder_tool_layout = QVBoxLayout()
        binder_tool_layout.setSpacing(5)
        binder_tool_layout.addWidget(self.netmhcpan_checkbox)
        binder_tool_layout.addWidget(self.mhcflurry_checkbox)
        binder_tool_layout.addWidget(self.bigmhc_checkbox)
        binder_tool_layout.addWidget(self.netmhcIIpan_checkbox)
        binder_tool_layout.addWidget(self.mixmhc2pred_checkbox)
        binder_layout.addLayout(binder_tool_layout)

        binder_rank_layout = QVBoxLayout()
        binder_rank_layout.setSpacing(5)
        netmhcpan_layout = QHBoxLayout()
        netmhcpan_layout.addWidget(netmhcpan_strong_label)
        netmhcpan_layout.addWidget(self.netmhcpan_strong_spinbox)
        netmhcpan_layout.addWidget(netmhcpan_weak_label)
        netmhcpan_layout.addWidget(self.netmhcpan_weak_spinbox)
        binder_rank_layout.addLayout(netmhcpan_layout)

        mhcflurry_layout = QHBoxLayout()
        mhcflurry_layout.addWidget(mhcflurry_strong_label)
        mhcflurry_layout.addWidget(self.mhcflurry_strong_spinbox)
        mhcflurry_layout.addWidget(mhcflurry_weak_label)
        mhcflurry_layout.addWidget(self.mhcflurry_weak_spinbox)
        binder_rank_layout.addLayout(mhcflurry_layout)

        bigmhc_layout = QHBoxLayout()
        bigmhc_layout.addWidget(bigmhc_strong_label)
        bigmhc_layout.addWidget(self.bigmhc_strong_spinbox)
        bigmhc_layout.addWidget(bigmhc_weak_label)
        bigmhc_layout.addWidget(self.bigmhc_weak_spinbox)
        binder_rank_layout.addLayout(bigmhc_layout)

        netmhcIIpan_layout = QHBoxLayout()
        netmhcIIpan_layout.addWidget(netmhcIIpan_strong_label)
        netmhcIIpan_layout.addWidget(self.netmhcIIpan_strong_spinbox)
        netmhcIIpan_layout.addWidget(netmhcIIpan_weak_label)
        netmhcIIpan_layout.addWidget(self.netmhcIIpan_weak_spinbox)
        binder_rank_layout.addLayout(netmhcIIpan_layout)

        mixmhc2pred_layout = QHBoxLayout()
        mixmhc2pred_layout.addWidget(mixmhc2pred_strong_label)
        mixmhc2pred_layout.addWidget(self.mixmhc2pred_strong_spinbox)
        mixmhc2pred_layout.addWidget(mixmhc2pred_weak_label)
        mixmhc2pred_layout.addWidget(self.mixmhc2pred_weak_spinbox)
        binder_rank_layout.addLayout(mixmhc2pred_layout)

        binder_rank_widget = QWidget()
        binder_rank_widget.setFixedWidth(600)
        binder_rank_widget.setLayout(binder_rank_layout)
        binder_layout.addWidget(binder_rank_widget)
        binder_layout.addSpacing(100)
        pept_info_layout.addLayout(binder_layout)

        # Immunogenicity prediction
        self.ig_checkbox = QCheckBox('Immunogenicity prediction')
        ig_layout = QHBoxLayout()
        ig_layout.setAlignment(Qt.AlignLeft)
        ig_layout.addWidget(self.ig_checkbox)
        pept_info_layout.addLayout(ig_layout)
        # Logo analysis

        pept_info_groupbox.setLayout(pept_info_layout)
        reporter_layout.addWidget(pept_info_groupbox)


        ### Protein
        prot_info_groupbox = QGroupBox('Protein Info')
        prot_info_layout = QVBoxLayout()
        prot_info_layout.insertSpacing(0, 5)
        prot_info_layout.setSpacing(20)

        protein_prophet_checkbox = QCheckBox('Protein inference (ProteinProphet)')
        protein_prophet_label = QLabel('CLI params:')
        self.protein_prophet_inputbox = QLineEdit()
        self.protein_prophet_inputbox.setText('--maxppmdiff 2000000')
        protein_prophet_layout = QHBoxLayout()
        protein_prophet_layout.setSpacing(5)
        protein_prophet_layout.addWidget(protein_prophet_label)
        protein_prophet_layout.addWidget(self.protein_prophet_inputbox)
        prot_infer_layout = QHBoxLayout()
        prot_infer_layout.setSpacing(50)
        prot_infer_layout.addWidget(protein_prophet_checkbox)
        prot_infer_layout.addLayout(protein_prophet_layout)
        prot_info_layout.addLayout(prot_infer_layout)

        fasta_checkbox = QCheckBox('Extract protein details from fasta')
        report_fasta_label = QLabel('FASTA path: ')
        self.report_fasta_inputbox = QLineEdit()
        self.report_fasta_inputbox.setPlaceholderText('Select .fasta or .fasta.fas file...')
        self.report_fasta_button = QPushButton("Select")
        self.report_fasta_button.clicked.connect(self.open_file_dialog)
        fasta_input_layout = QHBoxLayout()
        fasta_input_layout.setSpacing(5)
        fasta_input_layout.addWidget(report_fasta_label)
        fasta_input_layout.addWidget(self.report_fasta_inputbox)
        fasta_input_layout.addWidget(self.report_fasta_button)
        report_fasta_layout = QHBoxLayout()
        report_fasta_layout.setSpacing(50)
        report_fasta_layout.addWidget(fasta_checkbox)
        report_fasta_layout.addLayout(fasta_input_layout)
        prot_info_layout.addLayout(report_fasta_layout)

        prot_info_groupbox.setLayout(prot_info_layout)
        reporter_layout.addWidget(prot_info_groupbox)


        ### Report
        report_groupbox = QGroupBox('Report')
        report_layout = QVBoxLayout()
        report_layout.insertSpacing(0, 5)
        report_layout.setSpacing(20)

        # FDR
        self.fdr_checkbox = QCheckBox('FDR filter')
        psm_fdr_label = QLabel('PSM FDR:')
        self.psm_fdr_spinbox = QDoubleSpinBox()
        self.psm_fdr_spinbox.setRange(0.000001, 1)
        self.psm_fdr_spinbox.setValue(0.01)
        self.psm_fdr_spinbox.setSingleStep(0.01)
        pep_fdr_label = QLabel('Peptide FDR:')
        self.pep_fdr_spinbox = QDoubleSpinBox()
        self.pep_fdr_spinbox.setRange(0.000001, 1)
        self.pep_fdr_spinbox.setValue(0.01)
        self.pep_fdr_spinbox.setSingleStep(0.01)
        prot_fdr_label = QLabel('Protein FDR:')
        self.prot_fdr_spinbox = QDoubleSpinBox()
        self.prot_fdr_spinbox.setRange(0.000001, 1)
        self.prot_fdr_spinbox.setValue(0.01)
        self.prot_fdr_spinbox.setSingleStep(0.01)

        fdr_layout = QHBoxLayout()
        fdr_layout.setSpacing(50)
        psm_fdr_layout = QHBoxLayout()
        psm_fdr_layout.setSpacing(5)
        psm_fdr_layout.addWidget(psm_fdr_label)
        psm_fdr_layout.addWidget(self.psm_fdr_spinbox)
        pep_fdr_layout = QHBoxLayout()
        pep_fdr_layout.setSpacing(5)
        pep_fdr_layout.addWidget(pep_fdr_label)
        pep_fdr_layout.addWidget(self.pep_fdr_spinbox)
        prot_fdr_layout = QHBoxLayout()
        prot_fdr_layout.setSpacing(5)
        prot_fdr_layout.addWidget(prot_fdr_label)
        prot_fdr_layout.addWidget(self.prot_fdr_spinbox)
        fdr_layout.addWidget(self.fdr_checkbox)
        fdr_layout.addLayout(psm_fdr_layout)
        fdr_layout.addLayout(pep_fdr_layout)
        fdr_layout.addLayout(prot_fdr_layout)
        fdr_layout.addSpacing(300)
        report_layout.addLayout(fdr_layout)
        report_groupbox.setLayout(report_layout)
        reporter_layout.addWidget(report_groupbox)

        # Contaminants
        self.contam_checkbox = QCheckBox('Remove contaminants')
        contam_fasta_label = QLabel('Custom contaminant FASTA (optional): ')
        self.contam_fasta_inputbox = QLineEdit()
        self.contam_fasta_inputbox.setPlaceholderText('Select .fasta or .fasta.fas file for custom contaminants...')
        self.contam_fasta_button = QPushButton("Select")
        self.contam_fasta_button.clicked.connect(self.open_file_dialog)
        contam_layout = QHBoxLayout()
        contam_layout.setSpacing(50)
        contam_fasta_layout = QHBoxLayout()
        contam_fasta_layout.setSpacing(5)
        contam_fasta_layout.addWidget(contam_fasta_label)
        contam_fasta_layout.addWidget(self.contam_fasta_inputbox)
        contam_fasta_layout.addWidget(self.contam_fasta_button)
        contam_layout.addWidget(self.contam_checkbox)
        contam_layout.addLayout(contam_fasta_layout)
        report_layout.addLayout(contam_layout)

        ### Logger
        self.reporter_log_output = QTextEdit()
        self.reporter_log_output.setReadOnly(True)
        # self.reporter_log_output.setFixedHeight(150)
        reporter_layout.addWidget(self.reporter_log_output)

        ### Execution
        reporter_open_button = QPushButton("Open result folder")
        reporter_run_button = QPushButton('Generate report')
        reporter_run_layout = QHBoxLayout()
        # reporter_run_layout.setAlignment(Qt.AlignRight)
        reporter_run_layout.addWidget(reporter_open_button)
        reporter_run_layout.addWidget(reporter_run_button)
        reporter_layout.addLayout(reporter_run_layout)

        reporter_tab = QWidget()
        reporter_tab.setLayout(reporter_layout)
        self.tab_widget.addTab(reporter_tab, 'Report')


    def open_folder_dialog(self):
        sender = self.sender()
        file_dialog = QFileDialog(self)
        file_dialog.setFileMode(QFileDialog.Directory)

        if file_dialog.exec():
            selected_path = file_dialog.selectedFiles()[0]
            if sender == self.psm_button:
                self.psm_inputbox.setText(selected_path)
            elif sender == self.mzml_button:
                self.mzml_inputbox.setText(selected_path)
            elif sender == self.psm_output_button:
                self.psm_output_inputbox.setText(selected_path)
                self.result_inputbox.setText(selected_path)
            elif sender == self.raw_button:
                self.raw_inputbox.setText(selected_path)
            elif sender == self.raw_output_button:
                self.raw_output_inputbox.setText(selected_path)
                self.result_inputbox.setText(selected_path)
            elif sender == self.result_button:
                self.result_inputbox.setText(selected_path)

    def open_file_dialog(self):
        sender = self.sender()
        file_dialog = QFileDialog(self)
        file_dialog.setFileMode(QFileDialog.ExistingFiles)

        if file_dialog.exec():
            selected_path = file_dialog.selectedFiles()[0]
            if sender == self.allele_button:
                self.allele_inputbox.setText(selected_path)
            elif sender == self.fasta_button:
                self.fasta_inputbox.setText(selected_path)
                self.report_fasta_inputbox.setText(selected_path)
            elif sender == self.report_fasta_button:
                self.report_fasta_inputbox.setText(selected_path)
            elif sender == self.contam_fasta_button:
                self.contam_fasta_inputbox.setText(selected_path)
            elif sender == self.param_button:
                self.param_inputbox.setText(selected_path)
            elif sender == self.msfragger_browse_button:
                self.msfragger_inputbox.setText(selected_path)
            elif sender == self.autort_browse_button:
                self.autort_inputbox.setText(selected_path)
            elif sender == self.bigmhc_browse_button:
                self.bigmhc_inputbox.setText(selected_path)
            elif sender == self.netmhcpan_browse_button:
                self.netmhcpan_inputbox.setText(selected_path)
            elif sender == self.netmhcIIpan_browse_button:
                self.netmhcIIpan_inputbox.setText(selected_path)
            elif sender == self.mixmhc2pred_browse_button:
                self.mixmhc2pred_inputbox.setText(selected_path)


    def on_mhc_I_checkbox_toggled(self, checked):
        if checked:
            for checkbox in self.checkboxes_mhc_II:
                checkbox.setChecked(False)

    def on_mhc_II_checkbox_toggled(self, checked):
        if checked:
            for checkbox in self.checkboxes_mhc_I:
                checkbox.setChecked(False)

    def on_autopred_checkbox_toggled(self, checked):
        if checked:
            for checkbox in self.checkboxes_rt + self.checkboxes_ms2 + self.checkboxes_ccs + self.checkboxes_pe:
                checkbox.setDisabled(True)
        else:
            for checkbox in self.checkboxes_rt + self.checkboxes_ms2 + self.checkboxes_ccs + self.checkboxes_pe:
                checkbox.setDisabled(False)

    def show_message(self, message):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Information)
        msg.setWindowTitle("Information")
        msg.setText(message)
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        msg.exec()

    def load_default_config(self):
        pass

    def save_config(self):
        pass

    def load_config(self):
        pass

    def save_params(self):
        date = datetime.now().strftime("%y_%m_%d")
        param_filename = f'mhcbooster-{date}.params'
        with open(Path(self.output_inputbox.text()) / param_filename, 'w') as f:
            f.write('')

    def on_exec_clicked(self):
        if self.button_run.text() == 'RUN':
            self.run()
        else:
            self.worker_stop()
            self.add_log('Process terminated.')

    def on_install_clicked(self):
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        if Path(self.msfragger_inputbox.text()).exists() and self.msfragger_inputbox.text().endswith('.zip'):
            install_msfragger(self.msfragger_inputbox.text())
        self.progress_bar.setValue(1)
        if Path(self.autort_inputbox.text()).exists() and self.autort_inputbox.text().endswith('.zip'):
            install_autort(self.autort_inputbox.text())
        self.progress_bar.setValue(2)
        if Path(self.bigmhc_inputbox.text()).exists() and self.bigmhc_inputbox.text().endswith('.zip'):
            install_bigmhc(self.bigmhc_inputbox.text())
        self.progress_bar.setValue(3)
        if Path(self.netmhcpan_inputbox.text()).exists() and self.netmhcpan_inputbox.text().endswith('.gz'):
            install_netmhcpan(self.netmhcpan_inputbox.text())
        self.progress_bar.setValue(4)
        if Path(self.netmhcIIpan_inputbox.text()).exists() and self.netmhcIIpan_inputbox.text().endswith('.gz'):
            install_netmhcIIpan(self.netmhcIIpan_inputbox.text())
        self.progress_bar.setValue(5)
        if Path(self.mixmhc2pred_inputbox.text()).exists() and self.mixmhc2pred_inputbox.text().endswith('.zip'):
            install_mixmhc2pred(self.mixmhc2pred_inputbox.text())
        self.progress_bar.setValue(6)
        self.refresh_third_party_status()
        self.progress_bar.setVisible(False)

    def refresh_third_party_status(self):
        third_party_folder = Path(__file__).parent.parent.parent/'third_party'
        for path in third_party_folder.iterdir():
            if 'AutoRT' in path.name:
                self.autort_inputbox.setText(str(path))
            if 'bigmhc' in path.name:
                self.bigmhc_inputbox.setText(str(path))
            if 'MixMHC2pred' in path.name:
                self.mixmhc2pred_inputbox.setText(str(path))
            if 'MSFragger' in path.name:
                self.msfragger_inputbox.setText(str(path))
            if 'netMHCIIpan' in path.name:
                self.netmhcIIpan_inputbox.setText(str(path))
            if 'netMHCpan' in path.name:
                self.netmhcpan_inputbox.setText(str(path))

    def worker_start(self):
        self.button_run.setText('STOP')
        self.worker_thread.start()

    def worker_stop(self):
        self.add_log(f'Stopping subprocess...')
        start_time = time.time()
        self.worker_thread.stop()
        print(time.time() - start_time)
        self.button_run.setText('RUN')


    def run(self):
        self.add_log(f'Running MhcBooster {__version__}...')
        # File
        psm_folder = Path(self.psm_inputbox.text())
        mzml_folder = self.mzml_inputbox.text()
        output_folder = Path(self.psm_output_inputbox.text()).resolve()
        output_folder.mkdir(parents=True, exist_ok=True)
        pin_files = list(psm_folder.rglob('*.pin'))
        if len(pin_files) == 0:
            self.show_message('No pin files found')
            return

        # Run params
        fine_tune = self.ft_checkbox.isChecked()
        auto_pred = self.ap_checkbox.isChecked()
        min_pep_length, max_pep_length = None, None
        if self.pl_checkbox.isChecked():
            min_pep_length = self.pl_min.value()
            max_pep_length = self.pl_max.value()
        koina_server_url = self.koina_inputbox.text()
        n_threads = self.spinbox_thread.value()

        # App score
        app_predictors = []
        mhc_class = None
        for checkbox in self.checkboxes_mhc_I:
            if checkbox.isChecked():
                app_predictors.append(checkbox.text())
                mhc_class = 'I'
        for checkbox in self.checkboxes_mhc_II:
            if checkbox.isChecked():
                app_predictors.append(checkbox.text())
                mhc_class = 'II'

        allele_param = self.allele_inputbox.text()
        if len(app_predictors) > 0 and len(allele_param) == 0:
            self.show_message('Input alleles cannot be empty')
            return

        # RT score
        rt_predictors = []
        for checkbox in self.checkboxes_rt:
            if checkbox.isChecked():
                rt_predictors.append(checkbox.text())

        # MS2 score
        ms2_predictors = []
        for checkbox in self.checkboxes_ms2:
            if isinstance(checkbox, QLabel):
                continue
            if checkbox.isChecked():
                ms2_predictors.append(checkbox.text())

        # CCS score
        ccs_predictors = []
        for checkbox in self.checkboxes_ccs:
            if checkbox.isChecked():
                ccs_predictors.append(checkbox.text())

        # PE
        pe = False
        for checkbox in self.checkboxes_pe:
            if checkbox.isChecked():
                pe = True

        app_predictor_param = ' '.join(app_predictors)
        rt_predictor_param = ' '.join(rt_predictors)
        ms2_predictor_param = ' '.join(ms2_predictors)
        ccs_predictor_param = ' '.join(ccs_predictors)
        cli_path = str(Path(__file__).parent/'mhcbooster_cli.py')
        command = f'python {cli_path} -n {n_threads}'
        if min_pep_length and max_pep_length:
            command += f' --min_pep_len {min_pep_length} --max_pep_len {max_pep_length}'
        if len(app_predictor_param) > 0 and len(allele_param) > 0:
            command += f' --app_predictors {app_predictor_param}'
            command += f' --alleles {allele_param}'
            command += f' --mhc_class {mhc_class}'
        if len(rt_predictor_param) > 0:
            command += f' --rt_predictors {rt_predictor_param}'
        if len(ms2_predictor_param) > 0:
            command += f' --ms2_predictors {ms2_predictor_param}'
        if len(ccs_predictor_param) > 0:
            command += f' --ccs_predictors {ccs_predictor_param}'
        if pe:
            command += f' --encode_peptide_sequences'
        if fine_tune:
            command += f' --fine_tune'
        if auto_pred:
            command += f' --auto_pred'
        if len(koina_server_url) > 0:
            command += f' --koina_server_url {koina_server_url}'
        command += f' --input {psm_folder} --output_dir {output_folder}'
        if len(mzml_folder) > 0:
            command += f' --mzml_dir {mzml_folder}'

        self.worker_thread.commands = [command]
        self.worker_start()


    def add_log(self, message):
        print(message)
        # if '\r' in message:
        #     print('Bingo!')
        #     self.log_output.moveCursor(QTextCursor.StartOfLine)
        self.log_output.append(message)
        self.log_output.moveCursor(QTextCursor.End)
        self.log_output.ensureCursorVisible()


class MhcBoosterWorker(QThread):
    message = Signal(str)
    finished = Signal()
    def __init__(self, commands):
        super().__init__()
        self.commands = commands
        self.process = None
        self._stop_flag = False

    def run(self):
        self._stop_flag = False
        for command in self.commands:
            if self._stop_flag:
                break

            self.process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            while True:
                if self._stop_flag:
                    self.process.terminate()  # Try to terminate the process gracefully
                    self.message.emit("Stopped gracefully...")
                    try:
                        self.process.wait(timeout=2) # Wait a bit to allow graceful termination
                    except subprocess.TimeoutExpired:
                        self.message.emit("Process didn't terminate in time, forcing kill...")
                        self.process.kill()
                    break

                    # Read output of the subprocess if needed
                stdout_line = self.process.stdout.readline()
                if stdout_line:
                    msg = stdout_line.decode('utf-8').strip()
                    self.message.emit(msg)

                if not stdout_line:
                    ret_code = self.process.poll()
                    if ret_code is not None:
                        self.message.emit(f"Process finished with return code: {ret_code}")
                        break
                time.sleep(0.01)  # Add a small delay to avoid high CPU usage
        self.finished.emit()

    def stop(self):
        self._stop_flag = True
        self.quit()  # Quit the QThread event loop
        # self.wait()  # Wait for the thread to finish

def run():
    app = QApplication(sys.argv)
    gui = MhcBoosterGUI()
    gui.show()

    sys.exit(app.exec())

if __name__ == '__main__':
    run()