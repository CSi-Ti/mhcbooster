from PySide2.QtGui import QTextCursor
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit
import sys
import time

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Carriage Return Handling")

        # Create a QTextEdit for displaying the log
        self.text_edit = QTextEdit()
        self.text_edit.setReadOnly(True)

        # Set up layout
        layout = QVBoxLayout()
        layout.addWidget(self.text_edit)
        self.setLayout(layout)

        # Example of updating with \r
        self.log_with_carriage_return()

    def log_with_carriage_return(self):
        """Method that adds text with carriage return to the QTextEdit."""
        self.text_edit.append('Start')
        for i in range(5):
            # Simulate output with \r (this will overwrite the line)
            text = f"Loading {i}...\r"
            if text.endswith("\r"):
                self.text_edit.moveCursor(QTextCursor.StartOfBlock)
            self.text_edit.append(text)  # Add text, it will overwrite in a loop
            time.sleep(1)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
