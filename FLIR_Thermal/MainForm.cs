﻿using System;
using System.Diagnostics;
using System.Windows.Forms;
using Flir.Atlas.Live.Device;
using Flir.Atlas.Live.Discovery;
using System.Drawing;
using System.Net;
using System.Net.Sockets;
using System.Text;

namespace DerekCam
{
    public partial class MainForm : Form
    {
        int gmrCount = 0;
        private Camera _camera;
        readonly Timer _timerRefreshUi = new Timer();
        public MainForm()
        {
            InitializeComponent();
            _timerRefreshUi.Interval = 20;
            _timerRefreshUi.Tick += _timerRefreshUi_Tick;
            StartClient();
        }

        public static void StartClient()
        {
            // Data buffer for incoming data.  
            byte[] bytes = new byte[1024];

            // Connect to a remote device.  
            try
            {
                // Establish the remote endpoint for the socket.  
                // This example uses port 8078 on the local computer.  
                //IPHostEntry ipHostInfo = Dns.GetHostEntry(Dns.GetHostName());
                IPAddress ipAddress = System.Net.IPAddress.Parse("192.168.1.127");
                IPEndPoint remoteEP = new IPEndPoint(ipAddress, 8078);

                // Create a TCP/IP  socket.  
                Socket sender = new Socket(ipAddress.AddressFamily,
                    SocketType.Stream, ProtocolType.Tcp);

                // Connect the socket to the remote endpoint. Catch any errors.  
                try
                {
                    sender.Connect(remoteEP);

                    Console.WriteLine("Socket connected to {0}",
                        sender.RemoteEndPoint.ToString());

                    // Encode the data string into a byte array.  
                    byte[] msg = Encoding.ASCII.GetBytes("This is a test<EOF>");

                    // Send the data through the socket.  
                    int bytesSent = sender.Send(msg);

                    // Receive the response from the remote device.  
                    int bytesRec = sender.Receive(bytes);
                    Console.WriteLine("Echoed test = {0}",
                        Encoding.ASCII.GetString(bytes, 0, bytesRec));

                    // Release the socket.  
                    sender.Shutdown(SocketShutdown.Both);
                    sender.Close();

                }
                catch (ArgumentNullException ane)
                {
                    Console.WriteLine("ArgumentNullException : {0}", ane.ToString());
                }
                catch (SocketException se)
                {
                    Console.WriteLine("SocketException : {0}", se.ToString());
                }
                catch (Exception e)
                {
                    Console.WriteLine("Unexpected exception : {0}", e.ToString());
                }

            }
            catch (Exception e)
            {
                Console.WriteLine(e.ToString());
            }
        }

        void _timerRefreshUi_Tick(object sender, EventArgs e)
        {
            if (!IsDirty) return;
            IsDirty = false;
            if (_camera == null) return;

            _camera.GetImage().EnterLock();
            try
            {
                pictureBox1.Image = _camera.GetImage().Image;
                //this is where we want to add in the machine learning
                gmrCount = gmrCount + 1;
                System.IO.File.WriteAllText(@"C:\labs\gmrTest.txt", gmrCount.ToString());
                if (gmrCount == 100) {
                    Bitmap img = new Bitmap(_camera.GetImage().Image);
                    img.Save(@"C:\labs\gmrImage3.bmp");
                }
            }
            catch (Exception exception)
            {
                Trace.TraceError(exception.Message);
            }
            finally
            {
                _camera.GetImage().ExitLock();
            }
        }

        private void discoveryToolStripMenuItem_Click(object sender, EventArgs e)
        {
            var discoveryDlg = new DiscoveryForm();
            if (discoveryDlg.ShowDialog() == DialogResult.OK)
            {
                ConnectCamera(discoveryDlg.SelectedCameraDevice);
            }
        }

        private void ConnectCamera(CameraDeviceInfo cameraDeviceInfo)
        {
            DisposeCamera();
            switch (cameraDeviceInfo.SelectedStreamingFormat)
            {
                case ImageFormat.FlirFileFormat:
                    _camera = new ThermalCamera();
                    break;
                case ImageFormat.Argb:
                    _camera = new VideoOverlayCamera();
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
            }
            _camera.ConnectionStatusChanged += _camera_ConnectionStatusChanged;
            _camera.GetImage().Changed += Image_Changed;
            _camera.Connect(cameraDeviceInfo);
            if (_camera.Recorder == null && _recorder != null)
            {
                _recorder.Dispose();
                _recorder = null;
            }
            if (_recorder != null)
            {
                if (!_recorder.IsDisposed)
                    _recorder.Initialize(_camera);
            }
            recorderToolStripMenuItem.Enabled = _camera.Recorder != null;
            _timerRefreshUi.Start();
        }

        private void DisposeCamera()
        {
            _timerRefreshUi.Stop();
            if (_camera == null) return;
            if (_recorder != null)
            {
                _recorder.UnInitialize();
                _recorder.Dispose();
            }
            _camera.ConnectionStatusChanged -= _camera_ConnectionStatusChanged;
            _camera.GetImage().Changed -= Image_Changed;
            _camera.Dispose();
        }

        private bool IsDirty { get; set; }
        void Image_Changed(object sender, Flir.Atlas.Image.ImageChangedEventArgs e)
        {
            IsDirty = true;
        }

        void _camera_ConnectionStatusChanged(object sender, Flir.Atlas.Live.ConnectionStatusChangedEventArgs e)
        {
            BeginInvoke((Action) (()=> toolStripConnectionStatus.Text = e.Status.ToString()));
            BeginInvoke((Action)(() => cameraToolStripMenuItem.Enabled = e.Status == ConnectionStatus.Connected)); 
        }

        private void DisconnectCamera()
        {
            if (_camera == null) return;
            
            _camera.Disconnect();
        }

        private void MainForm_Load(object sender, EventArgs e)
        {
            
        }

        private void disconnectToolStripMenuItem_Click(object sender, EventArgs e)
        {
            DisconnectCamera();
        }

        private void MainForm_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (_recorder != null && _recorder.IsDisposed == false)
            {
                _recorder.Dispose();    
            }
            
            DisposeCamera();
        }

        private RecorderForm _recorder;
        private void recorderToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (_recorder == null || _recorder.IsDisposed)
            {
                _recorder = new RecorderForm();
                _recorder.Initialize(_camera);
                _recorder.SelectedFileMouseDoubleClick += _recorder_SelectedFileMouseDoubleClick;
            }
            _recorder.Show();
            _recorder.Focus();
        }

        void _recorder_SelectedFileMouseDoubleClick(object sender, SelectedFileEventArgs e)
        {
            var playback = new PlaybackForm(e.FilePath);
            playback.Show();
        }

        private const string FileFilter = "IR Image Files(*.jpg;*.img;*.seq;*.fcf)|*.jpg;*.img;*.seq;*.fcf|All files (*.*)|*.*";

        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            using (var dialog = new OpenFileDialog())
            {
                dialog.InitialDirectory = Environment.GetFolderPath(Environment.SpecialFolder.MyPictures) + @"\FLIR";
                dialog.Filter = FileFilter;
                if (dialog.ShowDialog() != DialogResult.OK) return;
                var playback = new PlaybackForm(dialog.FileName);
                playback.Show();
            }
        }
    }
}
