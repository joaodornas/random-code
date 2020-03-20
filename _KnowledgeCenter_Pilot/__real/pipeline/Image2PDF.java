/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package image2pdf;

import java.io.FileInputStream;
import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.edit.PDPageContentStream;
import org.apache.pdfbox.pdmodel.graphics.xobject.PDJpeg;
import org.apache.pdfbox.pdmodel.graphics.xobject.PDXObjectImage;
import org.apache.pdfbox.pdmodel.common.*;

/**
 *
 * @author joaodornas
 */
public class Image2PDF {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        File folder = new File("/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/__data/DOCS_PAPERS/__encrypted/asIMAGE-PDFBOX-300/");
        File[] listOfFiles = folder.listFiles();

        for (int i = 0; i < listOfFiles.length; i++) {
        
        String imageName = listOfFiles[i].getName();
        String fileName = listOfFiles[i].getName() + ".pdf";
        
        try {

            PDDocument doc = new PDDocument();
            
            InputStream in = new FileInputStream(folder + "/" + imageName);
            BufferedImage bimg = ImageIO.read(in);
            float width = bimg.getWidth();
            float height = bimg.getHeight();
            PDPage page = new PDPage(new PDRectangle(width, height));
            
            doc.addPage(page);
            
            PDXObjectImage img = new PDJpeg(doc, new FileInputStream(folder + "/" + imageName));
            PDPageContentStream contentStream = new PDPageContentStream(doc, page);
            contentStream.drawImage(img, 0, 0);
            contentStream.close();
            in.close();

            doc.save(fileName);
            doc.close();

        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
            
    }

    }
}
        