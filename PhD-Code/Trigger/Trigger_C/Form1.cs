using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO.Ports;
using System.Net;
using System.Net.Sockets;
using System.IO;
using System.Diagnostics;

namespace Trigger
{
    public partial class triggerForm : Form
    {
        public string subjectName;

        public bool firstTR = true;
        public bool stopClicked = false;
        public DateTime startTime;

        public Socket sock = new Socket(AddressFamily.InterNetwork, SocketType.Dgram, ProtocolType.Udp);
        public IPAddress host = IPAddress.Parse("127.0.0.1");

        public Int32 udpPort;

        public StreamWriter fileText;
        public StreamWriter fileDate;

        public string fcnName;
        public string paramFcn;

        public triggerForm()
        {
            InitializeComponent();
        }

        private void triggerForm_Load(object sender, EventArgs e)
        {

            var serialPortSource = new List<comboBoxValues>();

            serialPortSource.Add(new comboBoxValues () { Name = "COM1", Value = "COM1" });
            serialPortSource.Add(new comboBoxValues () { Name = "COM2", Value = "COM2" });
            serialPortSource.Add(new comboBoxValues () { Name = "COM3", Value = "COM3" });
            serialPortSource.Add(new comboBoxValues () { Name = "COM4", Value = "COM4" });

            this.serialPortComboBox.DataSource = serialPortSource;
            this.serialPortComboBox.DisplayMember = "Name";
            this.serialPortComboBox.ValueMember = "Value";

            this.serialPortComboBox.DropDownStyle = ComboBoxStyle.DropDownList;

            var baudRateSource = new List<comboBoxValues>();

            baudRateSource.Add(new comboBoxValues() { Name = "9600", Value = "9600" });
            baudRateSource.Add(new comboBoxValues() { Name = "19200", Value = "19200" });
            baudRateSource.Add(new comboBoxValues() { Name = "38400", Value = "38400" });
            baudRateSource.Add(new comboBoxValues() { Name = "57600", Value = "57600" });
            baudRateSource.Add(new comboBoxValues() { Name = "115200", Value = "115200" });

            this.baudRateComboBox.DataSource = baudRateSource;
            this.baudRateComboBox.DisplayMember = "Name";
            this.baudRateComboBox.ValueMember = "Value";

            this.baudRateComboBox.DropDownStyle = ComboBoxStyle.DropDownList;

            this.udpPortTextBox.Text = "8080";

            this.fcnNameTextBox.Text = "Main";

        }

       public class comboBoxValues 
       {

           public string Name { get; set; }
           public string Value { get; set; }

       }

        private void startButton_Click(object sender, EventArgs e)
        {

            startTime = DateTime.Now;

            subjectName = this.subjectTextBox.Text;
  
            string comPort;
            comPort = this.serialPortComboBox.ValueMember;

            Int32 baudRate;
            baudRate = Convert.ToInt32(this.baudRateComboBox.ValueMember);

            Int16 dataBits;
            dataBits = Convert.ToInt16(this.dataBitsTextBox);  

            udpPort = Convert.ToInt32(this.udpPortTextBox);

            fcnName = this.fcnNameTextBox.Text;

            paramFcn = this.paramFcnTextBox.Text;

            SerialPort serialPort = new SerialPort();
            serialPort.PortName = comPort;
            serialPort.BaudRate = baudRate;
            serialPort.DataBits = dataBits;

            serialPort.Open();

            fileText = new StreamWriter(@"C:\INDIREA\" + subjectName + "-" + TRtime.Date + "-triggerDaatText.txt", true);
            fileDate = new StreamWriter(@"C:\INDIREA\" + subjectName + "-" + TRtime.Date + "-triggerDataDate.txt", true);

            while(true)
            {
                if (stopClicked);
                {
                    serialPort.Close();
                    break;
                }

                serialPort.PinChanged += new SerialPinChangedEventHandler(captureTrigger); 

            }

            }

        public void captureTrigger(object sender, SerialPinChangedEventArgs e)
        {

            DateTime TRtime = DateTime.Now;

            if (firstTR)
            {

                ProcessStartInfo startInfo = new ProcessStartInfo("matlab -nodesktop -nosplash -r " + fcnName + "(" + paramFcn + ").m &");
                startInfo.WindowStyle = ProcessWindowStyle.Minimized;

                Process.Start(startInfo);

                firstTR = false;

            }

            IPEndPoint endPoint = new IPEndPoint(host, udpPort);

            string triggerText = "TR";
            byte[] triggerBuffer = Encoding.ASCII.GetBytes(triggerText);

            sock.SendTo(triggerBuffer, endPoint);

            string dateTimeText;
            dateTimeText = TRtime.ToString("o");

            fileText.WriteLine(dateTimeText);
            fileDate.WriteLine(TRtime);

        }

        private void stopButton_Click(object sender, EventArgs e)
        {
            stopClicked = true;
        }


        }
}
