import sys
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout, QGroupBox, QLineEdit, QCheckBox, QFormLayout

class MyWindow(QWidget):
    def __init__(self):
        super().__init__()

        # 设置窗口标题和大小
        self.setWindowTitle("包围框示例")
        self.setGeometry(100, 100, 400, 300)

        # 创建主布局
        layout = QVBoxLayout(self)

        # 创建第一个包围框：包含一个输入框
        groupBox1 = QGroupBox("输入框")
        groupBox1Layout = QFormLayout(groupBox1)
        inputLine = QLineEdit()
        inputLine.setPlaceholderText("请输入文本...")
        groupBox1Layout.addRow("文本输入:", inputLine)

        # 创建第二个包围框：包含多个复选框
        groupBox2 = QGroupBox("多选框")
        groupBox2Layout = QVBoxLayout(groupBox2)
        checkBox1 = QCheckBox("选项 1")
        checkBox2 = QCheckBox("选项 2")
        checkBox3 = QCheckBox("选项 3")
        groupBox2Layout.addWidget(checkBox1)
        groupBox2Layout.addWidget(checkBox2)
        groupBox2Layout.addWidget(checkBox3)

        # 将包围框添加到主布局中
        layout.addWidget(groupBox1)
        layout.addWidget(groupBox2)

        self.setLayout(layout)

# 创建应用并运行
app = QApplication(sys.argv)
window = MyWindow()
window.show()
sys.exit(app.exec_())
