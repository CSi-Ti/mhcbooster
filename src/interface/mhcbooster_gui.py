# pip install PySide2
import os
import re
import subprocess
import sys
import logging
import time
from datetime import datetime
from pathlib import Path

from PySide2.QtCore import Qt, QThread, Signal
from PySide2.QtGui import QPixmap, QIcon, QTextCursor
from PySide2.QtWidgets import (QApplication, QWidget, QPushButton, QLabel,
                               QVBoxLayout, QLineEdit, QFileDialog, QHBoxLayout,
                               QCheckBox, QGridLayout, QSpinBox, QGroupBox,
                               QMessageBox, QTextEdit)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mhcvalidator import __version__


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
                    border: 1px solid lightgray;
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
            """)

        # GUI window
        self.setWindowTitle('MhcBooster')
        self.setWindowIcon(QIcon(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'caronlab_icon.png')))
        # self.setGeometry(100, 100, 800, 600)
        layout = QVBoxLayout()
        layout.setContentsMargins(50, 20, 50, 10) # left, top, right, bottom
        layout.setSpacing(30)

        ### INTRODUCTION
        logo_lab_label = QLabel()
        logo_pix_map = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'caronlab.png')).scaled(200, 150, Qt.KeepAspectRatio)
        logo_lab_label.setPixmap(logo_pix_map)
        logo_lab_label.resize(logo_pix_map.size())
        intro_label = QLabel('The Introduction of MhcBooster should be here. GitHub. Tutorial. Cite. CaronLab.')
        intro_layout = QHBoxLayout()
        intro_layout.addWidget(logo_lab_label)
        intro_layout.addWidget(intro_label)
        layout.addLayout(intro_layout)

        ### FILE MANAGEMENT
        file_groupbox = QGroupBox('Input / Output')
        file_group_layout = QVBoxLayout()
        file_group_layout.insertSpacing(0, 5)
        file_group_layout.setSpacing(5)

        # Input folder
        psm_path_label = QLabel('PSM folder: \t')
        self.psm_inputbox = QLineEdit()
        self.psm_inputbox.setPlaceholderText("Select input folder containing .pin files from Comet or MSFragger...")
        self.psm_button = QPushButton("Select")
        self.psm_button.clicked.connect(self.open_folder_dialog)
        psm_layout = QHBoxLayout()
        psm_layout.addWidget(psm_path_label)
        psm_layout.addWidget(self.psm_inputbox)
        psm_layout.addWidget(self.psm_button)
        file_group_layout.addLayout(psm_layout)

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
        file_group_layout.addLayout(mzml_layout)

        # Output folder
        output_label = QLabel('Output folder: \t')
        self.output_inputbox = QLineEdit()
        self.output_inputbox.setPlaceholderText('Select output folder...')
        self.output_button = QPushButton("Select")
        self.output_button.clicked.connect(self.open_folder_dialog)
        output_layout = QHBoxLayout()
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_inputbox)
        output_layout.addWidget(self.output_button)
        file_group_layout.addLayout(output_layout)
        file_groupbox.setLayout(file_group_layout)
        layout.addWidget(file_groupbox)

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
        layout.addWidget(mhc_groupbox)

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
        unsuitable_ms2_models = ['UniSpec', 'ms2pip_TTOF5600', 'Prosit_2024_intensity_XL_NMS2', 'Prosit_2023_intensity_XL_CMS2',
                                 'Prosit_2023_intensity_XL_CMS3']
        ms2_models = ['AlphaPeptDeep_ms2_generic',
                      'ms2pip_HCD2021', 'ms2pip_Immuno_HCD',
                      'ms2pip_timsTOF2023', 'ms2pip_timsTOF2024', 'ms2pip_iTRAQphospho', 'ms2pip_CID_TMT',
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
        layout.addWidget(gs_groupbox)


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
        layout.addWidget(rp_groupbox)

        ### Logger
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setFixedHeight(150)
        layout.addWidget(self.log_output)

        ### Execution
        self.button_run = QPushButton("RUN")
        self.button_run.clicked.connect(self.on_exec_clicked)
        layout.addWidget(self.button_run)

        self.worker_thread = MhcBoosterWorker(commands=None)
        self.worker_thread.message.connect(self.add_log)
        self.worker_thread.finished.connect(self.worker_stop)

        ### Footnote
        foot_label = QLabel('CaronLab 2024')
        foot_label.setAlignment(Qt.AlignRight)
        layout.addWidget(foot_label)

        self.setLayout(layout)

    def open_folder_dialog(self):
        sender = self.sender()
        file_dialog = QFileDialog(self)
        file_dialog.setFileMode(QFileDialog.Directory)

        if file_dialog.exec_():
            selected_path = file_dialog.selectedFiles()[0]
            if sender == self.psm_button:
                self.psm_inputbox.setText(selected_path)
            elif sender == self.mzml_button:
                self.mzml_inputbox.setText(selected_path)
            elif sender == self.output_button:
                self.output_inputbox.setText(selected_path)

    def open_file_dialog(self):
        sender = self.sender()
        file_dialog = QFileDialog(self)
        file_dialog.setFileMode(QFileDialog.ExistingFiles)

        if file_dialog.exec_():
            selected_path = file_dialog.selectedFiles()[0]
            if sender == self.allele_button:
                self.allele_inputbox.setText(selected_path)

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
        msg.exec_()

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
        pin_files = list(Path(self.psm_inputbox.text()).rglob('*.pin'))
        mzml_folder = self.mzml_inputbox.text()
        output_folder = Path(self.output_inputbox.text()).resolve()
        output_folder.mkdir(parents=True, exist_ok=True)
        if len(pin_files) == 0:
            self.show_message('No pin files found')

        # Run params
        fine_tune = self.ft_checkbox.isChecked()
        min_pep_length, max_pep_length = None, None
        if self.pl_checkbox.isChecked():
            min_pep_length = self.pl_min.value()
            max_pep_length = self.pl_max.value()
        koina_server_url = self.koina_inputbox.text()
        n_threads = self.spinbox_thread.value()

        # App score
        app_predictors = []
        for checkbox in self.checkboxes_mhc_I:
            if checkbox.isChecked():
                app_predictors.append(checkbox.text())
        for checkbox in self.checkboxes_mhc_II:
            if checkbox.isChecked():
                app_predictors.append(checkbox.text())

        alleles = []
        allele_map = {}
        if len(app_predictors) > 0:
            if self.allele_inputbox.text():
                if os.path.exists(self.allele_inputbox.text()):
                    for line in open(self.allele_inputbox.text()):
                        line_split = re.split(r'[\t,]', line)
                        allele_map[line_split[0].strip()] = [allele.strip() for allele in line_split[1].split(';')]
                else:
                    alleles = [allele.strip() for allele in self.allele_inputbox.text().split(';')]
            if len(alleles) == 0 and len(allele_map) == 0:
                self.show_message('Input alleles cannot be empty')

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

        commands = []
        for pin in pin_files:
            file_name = pin.stem
            print(file_name)
            run_alleles = alleles.copy()
            if len(alleles) == 0 and len(allele_map) != 0:
                for keyword in allele_map.keys():
                    if keyword in file_name:
                        run_alleles = allele_map[keyword]
                        break

            allele_param = ' '.join(run_alleles)
            app_predictor_param = ' '.join(app_predictors)
            rt_predictor_param = ' '.join(rt_predictors)
            ms2_predictor_param = ' '.join(ms2_predictors)
            ccs_predictor_param = ' '.join(ccs_predictors)
            command = f'python ../mhcvalidator/command_line.py -n {n_threads}'
            if min_pep_length and max_pep_length:
                command += f' --min_pep_len {min_pep_length} --max_pep_len {max_pep_length}'
            if len(app_predictor_param) > 0 and len(allele_param) > 0:
                command += f' --app_predictors {app_predictor_param}'
                command += f' --alleles {allele_param}'
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
            if len(koina_server_url) > 0:
                command += f' --koina_server_url {koina_server_url}'
            command += f' --input {pin} --output_dir {output_folder}'
            if len(mzml_folder) > 0:
                command += f' --mzml_dir {mzml_folder}'

            commands.append(command)
        self.worker_thread.commands = commands
        self.worker_start()


    def add_log(self, message):
        print(message)
        if '\r' in message:
            print('Bingo!')
            self.log_output.moveCursor(QTextCursor.StartOfLine)
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


if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = MhcBoosterGUI()

    gui.show()

    sys.exit(app.exec_())