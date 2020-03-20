namespace ReferencesPapers
{
    partial class MainForm
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
            this.components = new System.ComponentModel.Container();
            System.Windows.Forms.TreeNode treeNode1 = new System.Windows.Forms.TreeNode("Tree");
            this.treeView = new System.Windows.Forms.TreeView();
            this.contextMenuStripTreeView = new System.Windows.Forms.ContextMenuStrip(this.components);
            this.addNewCollectionToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.deleteCollectionToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.addNewToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.deleteToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.addNewToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.contextMenuStripTreeView.SuspendLayout();
            this.SuspendLayout();
            // 
            // treeView
            // 
            this.treeView.ContextMenuStrip = this.contextMenuStripTreeView;
            this.treeView.Location = new System.Drawing.Point(2, 90);
            this.treeView.Name = "treeView";
            treeNode1.Name = "Tree";
            treeNode1.Text = "Tree";
            this.treeView.Nodes.AddRange(new System.Windows.Forms.TreeNode[] {
            treeNode1});
            this.treeView.Size = new System.Drawing.Size(571, 812);
            this.treeView.TabIndex = 0;
            // 
            // contextMenuStripTreeView
            // 
            this.contextMenuStripTreeView.ImageScalingSize = new System.Drawing.Size(24, 24);
            this.contextMenuStripTreeView.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.addNewCollectionToolStripMenuItem,
            this.deleteCollectionToolStripMenuItem});
            this.contextMenuStripTreeView.Name = "contextMenuStripTreeView";
            this.contextMenuStripTreeView.Size = new System.Drawing.Size(241, 97);
            // 
            // addNewCollectionToolStripMenuItem
            // 
            this.addNewCollectionToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.addNewToolStripMenuItem,
            this.deleteToolStripMenuItem});
            this.addNewCollectionToolStripMenuItem.Name = "addNewCollectionToolStripMenuItem";
            this.addNewCollectionToolStripMenuItem.Size = new System.Drawing.Size(240, 30);
            this.addNewCollectionToolStripMenuItem.Text = "Collection...";
            // 
            // deleteCollectionToolStripMenuItem
            // 
            this.deleteCollectionToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.addNewToolStripMenuItem1});
            this.deleteCollectionToolStripMenuItem.Name = "deleteCollectionToolStripMenuItem";
            this.deleteCollectionToolStripMenuItem.Size = new System.Drawing.Size(240, 30);
            this.deleteCollectionToolStripMenuItem.Text = "PDF";
            // 
            // addNewToolStripMenuItem
            // 
            this.addNewToolStripMenuItem.Name = "addNewToolStripMenuItem";
            this.addNewToolStripMenuItem.Size = new System.Drawing.Size(252, 30);
            this.addNewToolStripMenuItem.Text = "Add New";
            // 
            // deleteToolStripMenuItem
            // 
            this.deleteToolStripMenuItem.Name = "deleteToolStripMenuItem";
            this.deleteToolStripMenuItem.Size = new System.Drawing.Size(252, 30);
            this.deleteToolStripMenuItem.Text = "Delete";
            // 
            // addNewToolStripMenuItem1
            // 
            this.addNewToolStripMenuItem1.Name = "addNewToolStripMenuItem1";
            this.addNewToolStripMenuItem1.Size = new System.Drawing.Size(252, 30);
            this.addNewToolStripMenuItem1.Text = "Add New";
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(9F, 20F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(2182, 900);
            this.Controls.Add(this.treeView);
            this.Name = "MainForm";
            this.Text = "MainForm";
            this.Load += new System.EventHandler(this.mainFormLoad);
            this.contextMenuStripTreeView.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.TreeView treeView;
        private System.Windows.Forms.ContextMenuStrip contextMenuStripTreeView;
        private System.Windows.Forms.ToolStripMenuItem addNewCollectionToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem deleteCollectionToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem addNewToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem deleteToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem addNewToolStripMenuItem1;
    }
}

