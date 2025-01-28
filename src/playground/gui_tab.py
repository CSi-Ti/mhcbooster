import sys
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QApplication, QWidget, QTabWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton

class TabExample(QWidget):
    def __init__(self):
        super().__init__()

        # Set the window title and size
        self.setWindowTitle("PySide2 Tab Example")
        self.resize(400, 300)

        # Create the QTabWidget
        self.tab_widget = QTabWidget()

        # Create the first tab page
        tab1 = QWidget()
        tab1_layout = QVBoxLayout()
        tab1_layout.addWidget(QLabel("This is Tab 1"))
        tab1_layout.addWidget(QLineEdit())  # Add a QLineEdit widget
        tab1.setLayout(tab1_layout)

        # Create the second tab page
        tab2 = QWidget()
        tab2_layout = QVBoxLayout()
        tab2_layout.addWidget(QLabel("This is Tab 2"))
        tab2_layout.addWidget(QPushButton("Click me!"))  # Add a QPushButton widget
        tab2.setLayout(tab2_layout)

        # Add the tabs to the tab widget
        self.tab_widget.addTab(tab1, "Tab 1")
        self.tab_widget.addTab(tab2, "Tab 2")

        # Create the main layout and add the tab widget
        layout = QVBoxLayout()
        layout.addWidget(self.tab_widget)
        self.setLayout(layout)

# Create the application instance and window
app = QApplication(sys.argv)
window = TabExample()
window.show()

# Run the application
sys.exit(app.exec_())
