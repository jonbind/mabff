import sys
import os
import traceback
from datetime import datetime

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFormLayout, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QFileDialog, QPlainTextEdit, QProgressBar,
    QToolBar, QAction, QStatusBar, QMessageBox, QStyle,
)
from PyQt5.QtCore import Qt, QObject, pyqtSignal, QThread,QSize
from PyQt5.QtGui import QPalette, QColor

# External functions (assumed to be available)
from mAb_prepv3 import run_workflow


class GlycoPrepWorker(QObject):
    started = pyqtSignal()
    log = pyqtSignal(str)
    finished = pyqtSignal(bool, str)  # success, message

    def __init__(self, input_file: str, itp_output_file: str, molecule_name: str):
        super().__init__()
        self.input_file = input_file
        self.itp_output_file = itp_output_file
        self.molecule_name = molecule_name

    def run(self):
        self.started.emit()
        try:
            self.log.emit(f"Starting processing for: {self.input_file}")
            run_workflow(self.input_file, self.molecule_name, self.itp_output_file)
            self.log.emit("Processing complete.")
            self.finished.emit(True, "Done")
        except Exception as e:
            self.log.emit("An error occurred:\n" + traceback.format_exc())
            self.finished.emit(False, str(e))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GlycoPrep GUI")
        self.setMinimumSize(700, 500)
        self._thread = None
        self._worker = None

        # Central widget and main layout
        central = QWidget()
        self.setCentralWidget(central)
        self.main_layout = QVBoxLayout(central)
        self.main_layout.setContentsMargins(16, 16, 16, 16)
        self.main_layout.setSpacing(12)

        # Form section
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignRight)
        form.setFormAlignment(Qt.AlignTop)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)

        self.input_edit = QLineEdit()
        self.input_edit.setPlaceholderText("Select a .pdb input file...")
        self.input_edit.setClearButtonEnabled(True)
        self.input_button = QPushButton("Browse")
        self.input_button.clicked.connect(self.browse_input)

        input_row = QHBoxLayout()
        input_row.addWidget(self.input_edit, stretch=1)
        input_row.addWidget(self.input_button)
        form.addRow("Input File:", input_row)

        self.itp_edit = QLineEdit()
        self.itp_edit.setPlaceholderText("Select an output .itp file path...")
        self.itp_edit.setClearButtonEnabled(True)
        self.itp_button = QPushButton("Browse")
        self.itp_button.clicked.connect(self.browse_itp_file)

        itp_row = QHBoxLayout()
        itp_row.addWidget(self.itp_edit, stretch=1)
        itp_row.addWidget(self.itp_button)
        form.addRow("ITP Output File:", itp_row)

        self.molecule_name_edit = QLineEdit()
        self.molecule_name_edit.setPlaceholderText("Enter molecule name (e.g., GLYCO_01)")
        self.molecule_name_edit.setClearButtonEnabled(True)
        form.addRow("Molecule Name:", self.molecule_name_edit)

        # Run/Reset buttons
        btn_row = QHBoxLayout()
        self.run_button = QPushButton("Run")
        self.run_button.setIcon(self.style().standardIcon(QStyle.SP_MediaPlay))
        self.run_button.clicked.connect(self.run_script)

        self.reset_button = QPushButton("Reset")
        self.reset_button.setIcon(self.style().standardIcon(QStyle.SP_DialogResetButton))
        self.reset_button.clicked.connect(self.reset_fields)

        btn_row.addStretch(1)
        btn_row.addWidget(self.reset_button)
        btn_row.addWidget(self.run_button)

        # Progress bar
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        self.progress.setVisible(False)

        # Log window
        self.log_window = QPlainTextEdit()
        self.log_window.setReadOnly(True)
        self.log_window.setPlaceholderText("Logs will appear here...")
        self.log_window.setMinimumHeight(220)

        # Assemble layouts
        self.main_layout.addLayout(form)
        self.main_layout.addLayout(btn_row)
        self.main_layout.addWidget(self.progress)
        self.main_layout.addWidget(QLabel("Log:"))
        self.main_layout.addWidget(self.log_window, stretch=1)

        # Toolbar
        self._init_toolbar()

        # Status bar
        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self.status.showMessage("Ready")

        # Styling
        self._apply_modern_styles()

    def _init_toolbar(self):
        toolbar = QToolBar("Main Toolbar")
        toolbar.setIconSize(QSize(20, 20))
        self.addToolBar(toolbar)

        open_act = QAction(self.style().standardIcon(QStyle.SP_DialogOpenButton), "Open Input", self)
        open_act.triggered.connect(self.browse_input)
        toolbar.addAction(open_act)

        save_act = QAction(self.style().standardIcon(QStyle.SP_DialogSaveButton), "Select ITP Output", self)
        save_act.triggered.connect(self.browse_itp_file)
        toolbar.addAction(save_act)

        toolbar.addSeparator()

        run_act = QAction(self.style().standardIcon(QStyle.SP_MediaPlay), "Run", self)
        run_act.triggered.connect(self.run_script)
        toolbar.addAction(run_act)

        reset_act = QAction(self.style().standardIcon(QStyle.SP_DialogResetButton), "Reset", self)
        reset_act.triggered.connect(self.reset_fields)
        toolbar.addAction(reset_act)

    def _apply_modern_styles(self):
        # Rounded inputs and subtle hover for buttons
        self.setStyleSheet("""
            QLineEdit {
                padding: 6px 8px;
                border: 1px solid #666;
                border-radius: 6px;
                background: #2b2b2b;
                color: #f0f0f0;
            }
            QLineEdit:focus {
                border: 1px solid #4da3ff;
                background: #303030;
            }
            QPushButton {
                padding: 6px 14px;
                border: 1px solid #555;
                border-radius: 6px;
                background: #3a3a3a;
                color: #f0f0f0;
            }
            QPushButton:hover {
                background: #4a4a4a;
            }
            QPushButton:disabled {
                color: #999;
                background: #3a3a3a;
                border-color: #444;
            }
            QPlainTextEdit {
                border: 1px solid #666;
                border-radius: 6px;
                background: #222;
                color: #e6e6e6;
                font-family: Consolas, Menlo, Monaco, monospace;
                font-size: 12px;
            }
            QProgressBar {
                border: 1px solid #666;
                border-radius: 6px;
                background: #2b2b2b;
            }
            QProgressBar::chunk {
                background-color: #4da3ff;
                border-radius: 6px;
            }
        """)

    def browse_input(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select Input File", "", "PDB Files (*.pdb);;All Files (*)"
        )
        if file_name:
            self.input_edit.setText(file_name)
            self.status.showMessage(f"Selected input: {file_name}", 3000)

    def browse_itp_file(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Select ITP Output File", "", "ITP Files (*.itp);;All Files (*)"
        )
        if file_path:
            self.itp_edit.setText(file_path)
            self.status.showMessage(f"Selected output: {file_path}", 3000)

    def reset_fields(self):
        self.input_edit.clear()
        self.itp_edit.clear()
        self.molecule_name_edit.clear()
        self.log_window.clear()
        self.status.showMessage("Fields reset", 2000)

    def run_script(self):
        itp_output_file_path = self.itp_edit.text().strip()
        input_file_name = self.input_edit.text().strip()
        molecule_name = self.molecule_name_edit.text().strip()

        # Validate inputs
        errors = []
        if not input_file_name:
            errors.append("Please select an input file.")
        elif not os.path.isfile(input_file_name):
            errors.append("Input file does not exist.")
        elif not input_file_name.lower().endswith(".pdb"):
            errors.append("Input file must be a .pdb file.")

        if not itp_output_file_path:
            errors.append("Please select an ITP output file path.")
        elif not itp_output_file_path.lower().endswith(".itp"):
            errors.append("Output file should have a .itp extension.")

        if not molecule_name:
            errors.append("Please enter a molecule name.")

        if errors:
            self._show_error("\n".join(errors))
            return

        # Prepare logging and UI feedback
        self._append_log(f"Running with input: {input_file_name}")
        self._append_log(f"Output: {itp_output_file_path}")
        self._append_log(f"Molecule name: {molecule_name}")

        self._set_busy(True)

        # Create thread & worker
        self._thread = QThread()
        self._worker = GlycoPrepWorker(input_file_name, itp_output_file_path, molecule_name)
        self._worker.moveToThread(self._thread)

        # Connect signals
        self._thread.started.connect(self._worker.run)
        self._worker.started.connect(lambda: self.status.showMessage("Processing..."))
        self._worker.log.connect(self._append_log)
        self._worker.finished.connect(self._on_finished)
        self._thread.finished.connect(self._thread.deleteLater)

        # Start thread
        self._thread.start()

    def _append_log(self, msg: str):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_window.appendPlainText(f"[{timestamp}] {msg}")

    def _on_finished(self, success: bool, message: str):
        # Clean up thread
        if self._thread and self._thread.isRunning():
            self._thread.quit()
            self._thread.wait()
        self._thread = None
        self._worker = None

        self._set_busy(False)

        if success:
            self._append_log("Done!")
            self.status.showMessage("Completed successfully", 4000)
            QMessageBox.information(self, "Success", "Processing completed successfully.")
        else:
            self.status.showMessage("Error occurred", 4000)
            self._show_error(f"An error occurred:\n\n{message}")

    def _set_busy(self, busy: bool):
        self.run_button.setEnabled(not busy)
        self.progress.setVisible(busy)
        if busy:
            self.progress.setRange(0, 0)  # Indeterminate
            QApplication.setOverrideCursor(Qt.WaitCursor)
        else:
            self.progress.setRange(0, 100)
            QApplication.restoreOverrideCursor()

    def _show_error(self, text: str):
        QMessageBox.critical(self, "Error", text)
        self._append_log("ERROR: " + text)


def apply_dark_fusion(app: QApplication):
    app.setStyle("Fusion")
    dark_palette = QPalette()

    dark_color = QColor(45, 45, 45)
    disabled_color = QColor(127, 127, 127)
    text_color = QColor(220, 220, 220)
    highlight_color = QColor(77, 163, 255)

    dark_palette.setColor(QPalette.Window, dark_color)
    dark_palette.setColor(QPalette.WindowText, text_color)
    dark_palette.setColor(QPalette.Base, QColor(30, 30, 30))
    dark_palette.setColor(QPalette.AlternateBase, dark_color)
    dark_palette.setColor(QPalette.ToolTipBase, text_color)
    dark_palette.setColor(QPalette.ToolTipText, text_color)
    dark_palette.setColor(QPalette.Text, text_color)
    dark_palette.setColor(QPalette.Button, QColor(50, 50, 50))
    dark_palette.setColor(QPalette.ButtonText, text_color)
    dark_palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
    dark_palette.setColor(QPalette.Highlight, highlight_color)
    dark_palette.setColor(QPalette.HighlightedText, QColor(0, 0, 0))

    # Disabled
    dark_palette.setColor(QPalette.Disabled, QPalette.Text, disabled_color)
    dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, disabled_color)

    app.setPalette(dark_palette)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    apply_dark_fusion(app)

    # Optional: better scaling on HiDPI
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
