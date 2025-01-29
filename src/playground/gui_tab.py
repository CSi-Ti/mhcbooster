from PySide2.QtWidgets import QApplication, QWidget, QTabWidget, QVBoxLayout, QPushButton, QTextEdit, QLabel, QFrame
from PySide2.QtCore import Qt

app = QApplication([])

window = QWidget()

# Create the main layout
main_layout = QVBoxLayout(window)

# Create a QTabWidget
tab_widget = QTabWidget()

# Create some content for each tab
tab1_content = QLabel("This is content for Tab 1")
tab2_content = QLabel("This is content for Tab 2")

# Add some tabs with buttons
tab_widget.addTab(tab1_content, "Tab 1")
tab_widget.addTab(tab2_content, "Tab 2")

# Remove the border of the tab bar (no box around tabs)
tab_widget.setStyleSheet("""
    QTabWidget::pane {
        border: none;  /* Remove the tab box */
    }

    /* Add a line under the tab buttons */
    QTabBar::top {
        border-bottom: 2px solid #aaa;  /* Line under the entire tab bar */
    }
""")

# Add the QTabWidget to the layout
main_layout.addWidget(tab_widget)
line = QFrame()
line.setFrameShape(QFrame.HLine)  # This creates a horizontal line
line.setFrameShadow(QFrame.Sunken)  # Optional: gives a sunken effect to the line
main_layout.addWidget(line)  # Add the line directly below the tab bar
# Set the layout for the window
window.setLayout(main_layout)
window.setWindowTitle("Tab Bar without Box")
window.setFixedSize(400, 300)

window.show()

app.exec_()
