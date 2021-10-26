namespace DerekCam
{
    partial class RecorderForm
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
            this.columnHeader1 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.buttonPause = new System.Windows.Forms.Button();
            this.labelOutputPath = new System.Windows.Forms.Label();
            this.buttonRec = new System.Windows.Forms.Button();
            this.buttonBrowse = new System.Windows.Forms.Button();
            this.textBoxNumFramesPreRec = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.checkBoxPreRecording = new System.Windows.Forms.CheckBox();
            this.groupBoxPreRecording = new System.Windows.Forms.GroupBox();
            this.textBoxTimeSpan = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.radioButtonRecSpeed = new System.Windows.Forms.RadioButton();
            this.radioButtonRecSpeedInterval = new System.Windows.Forms.RadioButton();
            this.listViewRecordings = new System.Windows.Forms.ListView();
            this.groupBoxRecSpeed = new System.Windows.Forms.GroupBox();
            this.labelElapsedTime = new System.Windows.Forms.Label();
            this.labelLostImages = new System.Windows.Forms.Label();
            this.labelFrameCounter = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.labelStatus = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.groupBoxPreRecording.SuspendLayout();
            this.groupBoxRecSpeed.SuspendLayout();
            this.SuspendLayout();
            // 
            // columnHeader1
            // 
            this.columnHeader1.Text = "Recordings";
            this.columnHeader1.Width = 71;
            // 
            // buttonPause
            // 
            this.buttonPause.Enabled = false;
            this.buttonPause.Location = new System.Drawing.Point(200, 300);
            this.buttonPause.Name = "buttonPause";
            this.buttonPause.Size = new System.Drawing.Size(48, 23);
            this.buttonPause.TabIndex = 36;
            this.buttonPause.Text = "Pause";
            this.buttonPause.UseVisualStyleBackColor = true;
            this.buttonPause.Click += new System.EventHandler(this.buttonPause_Click);
            // 
            // labelOutputPath
            // 
            this.labelOutputPath.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.labelOutputPath.Location = new System.Drawing.Point(19, 328);
            this.labelOutputPath.Name = "labelOutputPath";
            this.labelOutputPath.Size = new System.Drawing.Size(229, 19);
            this.labelOutputPath.TabIndex = 23;
            this.labelOutputPath.Text = "Path";
            // 
            // buttonRec
            // 
            this.buttonRec.Enabled = false;
            this.buttonRec.Location = new System.Drawing.Point(146, 300);
            this.buttonRec.Name = "buttonRec";
            this.buttonRec.Size = new System.Drawing.Size(48, 23);
            this.buttonRec.TabIndex = 35;
            this.buttonRec.Text = "Rec";
            this.buttonRec.UseVisualStyleBackColor = true;
            this.buttonRec.Click += new System.EventHandler(this.buttonRec_Click);
            // 
            // buttonBrowse
            // 
            this.buttonBrowse.Location = new System.Drawing.Point(19, 300);
            this.buttonBrowse.Name = "buttonBrowse";
            this.buttonBrowse.Size = new System.Drawing.Size(60, 23);
            this.buttonBrowse.TabIndex = 34;
            this.buttonBrowse.Text = "Browse...";
            this.buttonBrowse.UseVisualStyleBackColor = true;
            this.buttonBrowse.Click += new System.EventHandler(this.buttonBrowse_Click);
            // 
            // textBoxNumFramesPreRec
            // 
            this.textBoxNumFramesPreRec.Location = new System.Drawing.Point(67, 44);
            this.textBoxNumFramesPreRec.Name = "textBoxNumFramesPreRec";
            this.textBoxNumFramesPreRec.Size = new System.Drawing.Size(55, 20);
            this.textBoxNumFramesPreRec.TabIndex = 2;
            this.textBoxNumFramesPreRec.Text = "30";
            this.textBoxNumFramesPreRec.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(17, 47);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(44, 13);
            this.label2.TabIndex = 1;
            this.label2.Text = "Frames:";
            // 
            // checkBoxPreRecording
            // 
            this.checkBoxPreRecording.AutoSize = true;
            this.checkBoxPreRecording.Location = new System.Drawing.Point(17, 24);
            this.checkBoxPreRecording.Name = "checkBoxPreRecording";
            this.checkBoxPreRecording.Size = new System.Drawing.Size(59, 17);
            this.checkBoxPreRecording.TabIndex = 0;
            this.checkBoxPreRecording.Text = "Enable";
            this.checkBoxPreRecording.UseVisualStyleBackColor = true;
            this.checkBoxPreRecording.CheckStateChanged += new System.EventHandler(this.checkBoxPreRecording_CheckStateChanged);

            // 
            // groupBoxPreRecording
            // 
            this.groupBoxPreRecording.Controls.Add(this.textBoxNumFramesPreRec);
            this.groupBoxPreRecording.Controls.Add(this.label2);
            this.groupBoxPreRecording.Controls.Add(this.checkBoxPreRecording);
            this.groupBoxPreRecording.Enabled = false;
            this.groupBoxPreRecording.Location = new System.Drawing.Point(19, 222);
            this.groupBoxPreRecording.Name = "groupBoxPreRecording";
            this.groupBoxPreRecording.Size = new System.Drawing.Size(229, 72);
            this.groupBoxPreRecording.TabIndex = 33;
            this.groupBoxPreRecording.TabStop = false;
            this.groupBoxPreRecording.Text = "Pre-recording";
            // 
            // textBoxTimeSpan
            // 
            this.textBoxTimeSpan.Location = new System.Drawing.Point(43, 65);
            this.textBoxTimeSpan.Name = "textBoxTimeSpan";
            this.textBoxTimeSpan.Size = new System.Drawing.Size(62, 20);
            this.textBoxTimeSpan.TabIndex = 3;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(108, 68);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(113, 13);
            this.label1.TabIndex = 4;
            this.label1.Text = "Time span in seconds.";
            // 
            // radioButtonRecSpeed
            // 
            this.radioButtonRecSpeed.AutoSize = true;
            this.radioButtonRecSpeed.Location = new System.Drawing.Point(12, 19);
            this.radioButtonRecSpeed.Name = "radioButtonRecSpeed";
            this.radioButtonRecSpeed.Size = new System.Drawing.Size(61, 17);
            this.radioButtonRecSpeed.TabIndex = 2;
            this.radioButtonRecSpeed.TabStop = true;
            this.radioButtonRecSpeed.Text = "Highest";
            this.radioButtonRecSpeed.UseVisualStyleBackColor = true;
            this.radioButtonRecSpeed.CheckedChanged += new System.EventHandler(this.radioButtonRecSpeed_CheckedChanged);
            // 
            // radioButtonRecSpeedInterval
            // 
            this.radioButtonRecSpeedInterval.AutoSize = true;
            this.radioButtonRecSpeedInterval.Location = new System.Drawing.Point(12, 42);
            this.radioButtonRecSpeedInterval.Name = "radioButtonRecSpeedInterval";
            this.radioButtonRecSpeedInterval.Size = new System.Drawing.Size(60, 17);
            this.radioButtonRecSpeedInterval.TabIndex = 2;
            this.radioButtonRecSpeedInterval.TabStop = true;
            this.radioButtonRecSpeedInterval.Text = "Interval";
            this.radioButtonRecSpeedInterval.UseVisualStyleBackColor = true;
            // 
            // listViewRecordings
            // 
            this.listViewRecordings.Columns.AddRange(new System.Windows.Forms.ColumnHeader[] {
            this.columnHeader1});
            this.listViewRecordings.Location = new System.Drawing.Point(173, 10);
            this.listViewRecordings.Name = "listViewRecordings";
            this.listViewRecordings.Size = new System.Drawing.Size(75, 108);
            this.listViewRecordings.TabIndex = 37;
            this.listViewRecordings.UseCompatibleStateImageBehavior = false;
            this.listViewRecordings.View = System.Windows.Forms.View.Details;
            this.listViewRecordings.MouseDoubleClick += new System.Windows.Forms.MouseEventHandler(this.listViewRecordings_MouseDoubleClick);
            // 
            // groupBoxRecSpeed
            // 
            this.groupBoxRecSpeed.Controls.Add(this.textBoxTimeSpan);
            this.groupBoxRecSpeed.Controls.Add(this.label1);
            this.groupBoxRecSpeed.Controls.Add(this.radioButtonRecSpeed);
            this.groupBoxRecSpeed.Controls.Add(this.radioButtonRecSpeedInterval);
            this.groupBoxRecSpeed.Enabled = false;
            this.groupBoxRecSpeed.Location = new System.Drawing.Point(19, 124);
            this.groupBoxRecSpeed.Name = "groupBoxRecSpeed";
            this.groupBoxRecSpeed.Size = new System.Drawing.Size(229, 92);
            this.groupBoxRecSpeed.TabIndex = 32;
            this.groupBoxRecSpeed.TabStop = false;
            this.groupBoxRecSpeed.Text = "Recording speed";
            // 
            // labelElapsedTime
            // 
            this.labelElapsedTime.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.labelElapsedTime.Location = new System.Drawing.Point(84, 91);
            this.labelElapsedTime.Name = "labelElapsedTime";
            this.labelElapsedTime.Size = new System.Drawing.Size(83, 23);
            this.labelElapsedTime.TabIndex = 24;
            this.labelElapsedTime.Text = "0";
            this.labelElapsedTime.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // labelLostImages
            // 
            this.labelLostImages.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.labelLostImages.Location = new System.Drawing.Point(84, 64);
            this.labelLostImages.Name = "labelLostImages";
            this.labelLostImages.Size = new System.Drawing.Size(83, 23);
            this.labelLostImages.TabIndex = 25;
            this.labelLostImages.Text = "0";
            this.labelLostImages.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // labelFrameCounter
            // 
            this.labelFrameCounter.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.labelFrameCounter.Location = new System.Drawing.Point(84, 37);
            this.labelFrameCounter.Name = "labelFrameCounter";
            this.labelFrameCounter.Size = new System.Drawing.Size(83, 23);
            this.labelFrameCounter.TabIndex = 26;
            this.labelFrameCounter.Text = "0";
            this.labelFrameCounter.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(16, 96);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(70, 13);
            this.label9.TabIndex = 27;
            this.label9.Text = "Elapsed time:";
            // 
            // labelStatus
            // 
            this.labelStatus.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.labelStatus.Location = new System.Drawing.Point(84, 10);
            this.labelStatus.Name = "labelStatus";
            this.labelStatus.Size = new System.Drawing.Size(83, 23);
            this.labelStatus.TabIndex = 28;
            this.labelStatus.Text = "Stopped";
            this.labelStatus.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(16, 69);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(66, 13);
            this.label6.TabIndex = 29;
            this.label6.Text = "Lost images:";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(16, 42);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(44, 13);
            this.label5.TabIndex = 30;
            this.label5.Text = "Images:";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(16, 15);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(40, 13);
            this.label3.TabIndex = 31;
            this.label3.Text = "Status:";
            // 
            // RecorderForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(264, 358);
            this.Controls.Add(this.buttonPause);
            this.Controls.Add(this.labelOutputPath);
            this.Controls.Add(this.buttonRec);
            this.Controls.Add(this.buttonBrowse);
            this.Controls.Add(this.groupBoxPreRecording);
            this.Controls.Add(this.listViewRecordings);
            this.Controls.Add(this.groupBoxRecSpeed);
            this.Controls.Add(this.labelElapsedTime);
            this.Controls.Add(this.labelLostImages);
            this.Controls.Add(this.labelFrameCounter);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.labelStatus);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label3);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedToolWindow;
            this.Name = "RecorderForm";
            this.Text = "Recorder";
            this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.RecorderForm_FormClosing);
            this.Load += new System.EventHandler(this.RecorderForm_Load);
            this.groupBoxPreRecording.ResumeLayout(false);
            this.groupBoxPreRecording.PerformLayout();
            this.groupBoxRecSpeed.ResumeLayout(false);
            this.groupBoxRecSpeed.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ColumnHeader columnHeader1;
        private System.Windows.Forms.Button buttonPause;
        private System.Windows.Forms.Label labelOutputPath;
        private System.Windows.Forms.Button buttonRec;
        private System.Windows.Forms.Button buttonBrowse;
        private System.Windows.Forms.TextBox textBoxNumFramesPreRec;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.CheckBox checkBoxPreRecording;
        private System.Windows.Forms.GroupBox groupBoxPreRecording;
        private System.Windows.Forms.TextBox textBoxTimeSpan;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.RadioButton radioButtonRecSpeed;
        private System.Windows.Forms.RadioButton radioButtonRecSpeedInterval;
        private System.Windows.Forms.ListView listViewRecordings;
        private System.Windows.Forms.GroupBox groupBoxRecSpeed;
        private System.Windows.Forms.Label labelElapsedTime;
        private System.Windows.Forms.Label labelLostImages;
        private System.Windows.Forms.Label labelFrameCounter;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.Label labelStatus;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label3;

    }
}