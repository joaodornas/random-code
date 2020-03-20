using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Xml;
using System.IO;

using static ReferencesPapers.Globals;

namespace ReferencesPapers
{
    

    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();


        }

        private void mainFormLoad(object sender, EventArgs e)
        {
            string xml = File.ReadAllText(Globals.appFolder + @"\tree.xml");

            XmlDocument xmlDoc = new XmlDocument();
            xmlDoc.LoadXml(xml);

            XmlNodeList dataNodes = xmlDoc.FirstChild.ChildNodes;

            foreach (XmlNode node in dataNodes)
            {
                TreeNode t = new TreeNode();
                t.Text = node.Name;
                t.Name = node.Name;

                TreeNodeCollection parent = treeView.Nodes;

                TreeNode[] parentNode = parent.Find("Tree", false);

                parentNode[0].Nodes.Add(t);
                  
            }

            treeView.ExpandAll();

            //int x = 0;
        }

       
    }
}
