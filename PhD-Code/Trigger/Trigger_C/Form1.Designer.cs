namespace Trigger
{
    partial class triggerForm
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.serialPortComboBox = new System.Windows.Forms.ComboBox();
            this.serialPortLabel = new System.Windows.Forms.Label();
            this.udpPortLabel = new System.Windows.Forms.Label();
            this.udpPortTextBox = new System.Windows.Forms.TextBox();
            this.fcnNameLabel = new System.Windows.Forms.Label();
            this.fcnNameTextBox = new System.Windows.Forms.TextBox();
            this.paramFcnLabel = new System.Windows.Forms.Label();
            this.paramFcnTextBox = new System.Windows.Forms.TextBox();
            this.startButton = new System.Windows.Forms.Button();
            this.stopButton = new System.Windows.Forms.Button();
            this.baudRateComboBox = new System.Windows.Forms.ComboBox();
            this.baudRatelabel = new System.Windows.Forms.Label();
            this.dataBitsLabel = new System.Windows.Forms.Label();
            this.dataBitsTextBox = new System.Windows.Forms.TextBox();
            this.subjectLabel = new System.Windows.Forms.Label();
            this.subjectTextBox = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // serialPortComboBox
            // 
            this.serialPortComboBox.FormattingEnabled = true;
            this.serialPortComboBox.Location = new System.Drawing.Point(133, 65);
            this.serialPortComboBox.Name = "serialPortComboBox";
            this.serialPortComboBox.Size = new System.Drawing.Size(145, 21);
            this.serialPortComboBox.TabIndex = 0;
            // 
            // serialPortLabel
            // 
            this.serialPortLabel.AutoSize = true;
            this.serialPortLabel.Location = new System.Drawing.Point(45, 73);
            this.serialPortLabel.Name = "serialPortLabel";
            this.serialPortLabel.Size = new System.Drawing.Size(55, 13);
            this.serialPortLabel.TabIndex = 1;
            this.serialPortLabel.Text = "Serial Port";
            // 
            // udpPortLabel
            // 
            this.udpPortLabel.AutoSize = true;
            this.udpPortLabel.Location = new System.Drawing.Point(45, 166);
            this.udpPortLabel.Name = "udpPortLabel";
            this.udpPortLabel.Size = new System.Drawing.Size(52, 13);
            this.udpPortLabel.TabIndex = 2;
            this.udpPortLabel.Text = "UDP Port";
            // 
            // udpPortTextBox
            // 
            this.udpPortTextBox.Location = new System.Drawing.Point(133, 159);
            this.udpPortTextBox.Name = "udpPortTextBox";
            this.udpPortTextBox.Size = new System.Drawing.Size(65, 20);
            this.udpPortTextBox.TabIndex = 3;
            // 
            // fcnNameLabel
            // 
            this.fcnNameLabel.AutoSize = true;
            this.fcnNameLabel.Location = new System.Drawing.Point(48, 197);
            this.fcnNameLabel.Name = "fcnNameLabel";
            this.fcnNameLabel.Size = new System.Drawing.Size(79, 13);
            this.fcnNameLabel.TabIndex = 4;
            this.fcnNameLabel.Text = "Function Name";
            // 
            // fcnNameTextBox
            // 
            this.fcnNameTextBox.Location = new System.Drawing.Point(133, 194);
            this.fcnNameTextBox.Name = "fcnNameTextBox";
            this.fcnNameTextBox.Size = new System.Drawing.Size(145, 20);
            this.fcnNameTextBox.TabIndex = 5;
            // 
            // paramFcnLabel
            // 
            this.paramFcnLabel.AutoSize = true;
            this.paramFcnLabel.Location = new System.Drawing.Point(48, 227);
            this.paramFcnLabel.Name = "paramFcnLabel";
            this.paramFcnLabel.Size = new System.Drawing.Size(60, 13);
            this.paramFcnLabel.TabIndex = 6;
            this.paramFcnLabel.Text = "Parameters";
            // 
            // paramFcnTextBox
            // 
            this.paramFcnTextBox.Location = new System.Drawing.Point(133, 227);
            this.paramFcnTextBox.Name = "paramFcnTextBox";
            this.paramFcnTextBox.Size = new System.Drawing.Size(145, 20);
            this.paramFcnTextBox.TabIndex = 7;
            // 
            // startButton
            // 
            this.startButton.Location = new System.Drawing.Point(81, 278);
            this.startButton.Name = "startButton";
            this.startButton.Size = new System.Drawing.Size(72, 53);
            this.startButton.TabIndex = 8;
            this.startButton.Text = "Start";
            this.startButton.UseVisualStyleBackColor = true;
            this.startButton.Click += new System.EventHandler(this.startButton_Click);
            // 
            // stopButton
            // 
            this.stopButton.Location = new System.Drawing.Point(189, 278);
            this.stopButton.Name = "stopButton";
            this.stopButton.Size = new System.Drawing.Size(72, 53);
            this.stopButton.TabIndex = 9;
            this.stopButton.Text = "Stop";
            this.stopButton.UseVisualStyleBackColor = true;
            this.stopButton.Click += new System.EventHandler(this.stopButton_Click);
            // 
            // baudRateComboBox
            // 
            this.baudRateComboBox.FormattingEnabled = true;
            this.baudRateComboBox.Location = new System.Drawing.Point(133, 92);
            this.baudRateComboBox.Name = "baudRateComboBox";
            this.baudRateComboBox.Size = new System.Drawing.Size(145, 21);
            this.baudRateComboBox.TabIndex = 10;
            // 
            // baudRatelabel
            // 
            this.baudRatelabel.AutoSize = true;
            this.baudRatelabel.Location = new System.Drawing.Point(45, 100);
            this.baudRatelabel.Name = "baudRatelabel";
            this.baudRatelabel.Size = new System.Drawing.Size(55, 13);
            this.baudRatelabel.TabIndex = 11;
            this.baudRatelabel.Text = "BaudRate";
            // 
            // dataBitsLabel
            // 
            this.dataBitsLabel.AutoSize = true;
            this.dataBitsLabel.Location = new System.Drawing.Point(45, 130);
            this.dataBitsLabel.Name = "dataBitsLabel";
            this.dataBitsLabel.Size = new System.Drawing.Size(47, 13);
            this.dataBitsLabel.TabIndex = 12;
            this.dataBitsLabel.Text = "DataBits";
            // 
            // dataBitsTextBox
            // 
            this.dataBitsTextBox.Location = new System.Drawing.Point(133, 123);
            this.dataBitsTextBox.Name = "dataBitsTextBox";
            this.dataBitsTextBox.Size = new System.Drawing.Size(65, 20);
            this.dataBitsTextBox.TabIndex = 13;
            // 
            // subjectLabel
            // 
            this.subjectLabel.AutoSize = true;
            this.subjectLabel.Location = new System.Drawing.Point(45, 43);
            this.subjectLabel.Name = "subjectLabel";
            this.subjectLabel.Size = new System.Drawing.Size(43, 13);
            this.subjectLabel.TabIndex = 14;
            this.subjectLabel.Text = "Subject";
            // 
            // subjectTextBox
            // 
            this.subjectTextBox.Location = new System.Drawing.Point(133, 36);
            this.subjectTextBox.Name = "subjectTextBox";
            this.subjectTextBox.Size = new System.Drawing.Size(145, 20);
            this.subjectTextBox.TabIndex = 15;
            // 
            // triggerForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(330, 361);
            this.Controls.Add(this.subjectTextBox);
            this.Controls.Add(this.subjectLabel);
            this.Controls.Add(this.dataBitsTextBox);
            this.Controls.Add(this.dataBitsLabel);
            this.Controls.Add(this.baudRatelabel);
            this.Controls.Add(this.baudRateComboBox);
            this.Controls.Add(this.stopButton);
            this.Controls.Add(this.startButton);
            this.Controls.Add(this.paramFcnTextBox);
            this.Controls.Add(this.paramFcnLabel);
            this.Controls.Add(this.fcnNameTextBox);
            this.Controls.Add(this.fcnNameLabel);
            this.Controls.Add(this.udpPortTextBox);
            this.Controls.Add(this.udpPortLabel);
            this.Controls.Add(this.serialPortLabel);
            this.Controls.Add(this.serialPortComboBox);
            this.Name = "triggerForm";
            this.Text = "Trigger";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ComboBox serialPortComboBox;
        private System.Windows.Forms.Label serialPortLabel;
        private System.Windows.Forms.Label udpPortLabel;
        private System.Windows.Forms.TextBox udpPortTextBox;
        private System.Windows.Forms.Label fcnNameLabel;
        private System.Windows.Forms.TextBox fcnNameTextBox;
        private System.Windows.Forms.Label paramFcnLabel;
        private System.Windows.Forms.TextBox paramFcnTextBox;
        private System.Windows.Forms.Button startButton;
        private System.Windows.Forms.Button stopButton;
        private System.Windows.Forms.ComboBox baudRateComboBox;
        private System.Windows.Forms.Label baudRatelabel;
        private System.Windows.Forms.Label dataBitsLabel;
        private System.Windows.Forms.TextBox dataBitsTextBox;
        private System.Windows.Forms.Label subjectLabel;
        private System.Windows.Forms.TextBox subjectTextBox;
    }
}

