﻿using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using Flir.Atlas.Live;
using Flir.Atlas.Live.Device;
using Flir.Atlas.Live.Recorder;

namespace DerekCam
{
    public partial class RecorderForm : Form
    {
        private Camera _stream;
        private readonly Stopwatch _elapedTimeRecording = new Stopwatch();
        private readonly Timer _timer = new Timer();
        public RecorderForm()
        {
            InitializeComponent();
        }

        /// <summary>
        /// Un-initialize the recorder component.
        /// </summary>
        public void UnInitialize()
        {
            _timer.Stop();
            if (_stream != null)
                _stream.ConnectionStatusChanged -= _stream_ConnectionStatusChanged;
            _stream = null;
        }

        bool IsInitalizing { get; set; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="camera"></param>
        public void Initialize(Camera camera)
        {
            if (camera == null)
            {
                throw new InvalidOperationException("camera is null");
            }
            _stream = camera;
            _stream.ConnectionStatusChanged += _stream_ConnectionStatusChanged;
            labelOutputPath.Text = Environment.GetFolderPath(Environment.SpecialFolder.MyPictures) + @"\FLIR\";

            _timer.Interval = 110;
            _timer.Tick += _timer_Tick;
            _timer.Start();
            UpdateStatus();
        }

        void _stream_ConnectionStatusChanged(object sender, ConnectionStatusChangedEventArgs e)
        {
            BeginInvoke((Action) UpdateStatus);
        }
        
        void UpdateStatus()
        {
            IsInitalizing = true;
            if (_stream != null && _stream.ConnectionStatus == ConnectionStatus.Connected)
            {
                if (_stream.Recorder.IsTimeSpanEnabled)
                {
                    radioButtonRecSpeedInterval.Checked = true;
                    textBoxTimeSpan.Enabled = true;
                }
                else
                {
                    radioButtonRecSpeed.Checked = true;
                    textBoxTimeSpan.Enabled = false;
                }
                if (_stream.Recorder is ThermalImageRecorder)
                {
                    groupBoxPreRecording.Enabled = true;
                    checkBoxPreRecording.Checked = (_stream.Recorder as ThermalImageRecorder).IsPreRecordingEnabled;
                    textBoxNumFramesPreRec.Text =
                        (_stream.Recorder as ThermalImageRecorder).NumberOfFramesToPreRecord.ToString(
                            CultureInfo.InvariantCulture);
                }
                groupBoxRecSpeed.Enabled = true;
                buttonPause.Enabled = true;
                buttonRec.Enabled = true;
                textBoxTimeSpan.Text = _stream.Recorder.TimeSpan.TotalSeconds.ToString(CultureInfo.InvariantCulture);
            }
            else
            {
                groupBoxPreRecording.Enabled = false;
                groupBoxRecSpeed.Enabled = false;
                buttonPause.Enabled = false;
                buttonRec.Enabled = false;
            }
            IsInitalizing = false;
        }

        void _timer_Tick(object sender, EventArgs e)
        {
            RefreshData();
        }

        void RefreshData()
        {
            if (_stream.Recorder == null) return;
            
            labelStatus.Text = _stream.Recorder.Status.ToString();
            labelFrameCounter.Text = _stream.Recorder.FrameCount.ToString(CultureInfo.InvariantCulture);
            labelLostImages.Text = _stream.Recorder.LostImages.ToString(CultureInfo.InvariantCulture);
            labelElapsedTime.Text = _elapedTimeRecording.Elapsed.ToString();
        }
        private void buttonRec_Click(object sender, EventArgs e)
        {
            if (_stream.Recorder.Status == RecorderState.Stopped || _stream.Recorder.Status == RecorderState.PreRecording)
            {
                // start recording
                if (radioButtonRecSpeedInterval.Checked && _stream.Recorder.IsTimeSpanEnabled == false)
                {
                    int sec;
                    if (int.TryParse(textBoxTimeSpan.Text, out sec) && sec > 0)
                    {
                        _stream.Recorder.EnableTimeSpan(TimeSpan.FromSeconds(sec));
                    }
                    else
                    {
                        MessageBox.Show("Error parsing time span.");
                        return;
                    }
                }
                else if (radioButtonRecSpeed.Checked && _stream.Recorder.IsTimeSpanEnabled)
                {
                    _stream.Recorder.DisableTimeSpan();
                }
               
                buttonRec.Text = "Stop";
                _stream.Recorder.Start(GetNextFileName());
                _elapedTimeRecording.Reset();
                _elapedTimeRecording.Start();
            }
            else
            {
                _stream.Recorder.Stop();
                buttonRec.Text = "Rec";
                _elapedTimeRecording.Stop();
                var lv = new ListViewItem(_stream.Recorder.FileName) {Tag = _stream.Recorder.FullPath};
                listViewRecordings.Items.Add(lv);
            }
        }

        private string GetNextFileName()
        {
            string fileName;
            do
            {
                fileName = CreateFileName();
            } while (System.IO.File.Exists(fileName));
            return fileName;
        }

        private int _nextIndex = 1;
        private string CreateFileName()
        {
            var fileName = labelOutputPath.Text;
            fileName += string.Format("{0:0000}", _nextIndex++);
            return fileName + _stream.Recorder.Extension;
        }

        private void buttonPause_Click(object sender, EventArgs e)
        {
            switch (_stream.Recorder.Status)
            {
                case RecorderState.Paused:
                    _stream.Recorder.Continue();
                    buttonPause.Text = "Pause";
                    break;
                case RecorderState.Recording:
                    _stream.Recorder.Pause();
                    buttonPause.Text = "Continue";
                    break;
            }
        }

        private void buttonBrowse_Click(object sender, EventArgs e)
        {
            var dialog = new FolderBrowserDialog();
            if (dialog.ShowDialog() == DialogResult.OK)
            {
                labelOutputPath.Text = dialog.SelectedPath;
            }
        }

        private void checkBoxPreRecording_CheckStateChanged(object sender, EventArgs e)
        {
            if (!(_stream.Recorder is ThermalImageRecorder)) return;
            if (checkBoxPreRecording.Checked)
            {
                int frames;
                if (int.TryParse(textBoxNumFramesPreRec.Text, out frames) && frames > 0)
                {
                    (_stream.Recorder as ThermalImageRecorder).EnablePreRecording(frames);
                }
                else
                {
                    MessageBox.Show("Error parsing number of frames in pre-recording.");
                }
            }
            else
            {
                (_stream.Recorder as ThermalImageRecorder).DisablePreRecording();
            }
        }

        /// <summary>
        /// Event fired when the selected file is double clicked.
        /// </summary>
        public event EventHandler<SelectedFileEventArgs> SelectedFileMouseDoubleClick;
        private void listViewRecordings_MouseDoubleClick(object sender, MouseEventArgs e)
        {
            var items = listViewRecordings.SelectedItems;
            if (items.Count > 0)
            {
                var lv = items[0];
                var fullPath = lv.Tag as string;
                OnMouseDoubleClick(new SelectedFileEventArgs(fullPath));
            }
        }

        void OnMouseDoubleClick(SelectedFileEventArgs args)
        {
            if (SelectedFileMouseDoubleClick != null)
            {
                SelectedFileMouseDoubleClick(this, args);
            }
        }

        private void radioButtonRecSpeed_CheckedChanged(object sender, EventArgs e)
        {
            if (IsInitalizing) return;
            if (radioButtonRecSpeed.Checked)
            {
                _stream.Recorder.DisableTimeSpan();
                textBoxTimeSpan.Enabled = false;
            }
            else
            { // time span
                textBoxTimeSpan.Enabled = true;
            }
        }

        private void RecorderForm_FormClosing(object sender, FormClosingEventArgs e)
        {
            e.Cancel = true;
            Hide();
        }

        private void RecorderForm_Load(object sender, EventArgs e)
        {
        
        }
    }

    /// <summary>
    /// Event Args for Selected file event
    /// </summary>
    public class SelectedFileEventArgs : EventArgs
    {
        internal SelectedFileEventArgs(string path)
        {
            FilePath = path;
        }
        /// <summary>
        /// The full path to the file.
        /// </summary>
        public string FilePath { get; private set; }
    }
}
