from PySide2.QtGui import QTextCursor
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit, QPushButton
import sys


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Move Cursor to Line Start")

        # Create a QTextEdit for displaying the log
        self.text_edit = QTextEdit(self)
        self.text_edit.setPlainText("Line 1\nLine 2\nLine 3\nLine 4")  # Sample text

        # Create a button to trigger the cursor movement
        self.move_button = QPushButton("Move Cursor to Start of Line 2", self)
        self.move_button.clicked.connect(self.move_cursor_to_line_start)

        # Set up layout
        layout = QVBoxLayout()
        layout.addWidget(self.text_edit)
        layout.addWidget(self.move_button)
        self.setLayout(layout)

    def move_cursor_to_line_start(self):
        """Move the cursor to the start of the second line."""
        cursor = self.text_edit.textCursor()

        # Move the cursor to the start of the second line
        # cursor.movePosition(cursor.StartOfBlock)  # Move cursor to the start of the current block

        # Set the cursor position in the QTextEdit
        self.text_edit.moveCursor(QTextCursor.StartOfLine)
        # self.text_edit.setTextCursor(cursor)

        # Optional: Highlight the start of the line (for visual feedback)
        self.text_edit.ensureCursorVisible()
        # self.text_edit.setTextCursor(cursor)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
