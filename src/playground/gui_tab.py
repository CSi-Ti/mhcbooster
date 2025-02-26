import sys
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QProgressBar, QPushButton

class AnimatedProgressBar(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Animated Progress Bar with Color Effect")

        # Create a layout
        layout = QVBoxLayout()

        # Create the progress bar
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setRange(0, 100)  # Set the range of the progress bar
        self.progress_bar.setValue(0)  # Start from 0
        self.progress_bar.setTextVisible(False)  # Hide text inside the progress bar
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid gray;
                border-radius: 5px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
                border-radius: 3px;
            }
        """)
        layout.addWidget(self.progress_bar)

        # Button to start the progress animation
        self.start_button = QPushButton("Start Progress", self)
        self.start_button.clicked.connect(self.start_progress)
        layout.addWidget(self.start_button)

        # Set layout for the window
        self.setLayout(layout)

        # Set up a timer for animation
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_progress)

        # Initial progress value
        self.progress_value = 0

    def start_progress(self):
        # Start the timer to update the progress bar value
        self.timer.start(50)  # 50 milliseconds interval to create smooth animation

    def update_progress(self):
        # Update the progress value smoothly (simulate animation)
        self.progress_value += 1
        if self.progress_value <= 100:
            self.progress_bar.setValue(self.progress_value)
            self.progress_bar.setStyleSheet(f"""
                QProgressBar::chunk {{
                    background-color: rgb({self.progress_value*2}, {255 - self.progress_value*2}, 100);
                    border-radius: 3px;
                }}
            """)
        else:
            # Stop the animation when the progress reaches 100
            self.timer.stop()
            self.progress_bar.setStyleSheet("""
                QProgressBar {
                    border: 2px solid gray;
                    border-radius: 5px;
                    text-align: center;
                }
                QProgressBar::chunk {
                    background-color: #4CAF50;
                    border-radius: 3px;
                }
            """)
            self.start_button.setText("Animation Completed")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AnimatedProgressBar()
    window.show()
    sys.exit(app.exec())
