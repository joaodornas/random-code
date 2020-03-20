<?xml version="1.0"?>

<!-- xml2freemind.xsl $Id$
     Author: Terry Brown
     Created: Wed May 04 2005
-->

<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xi="http://www.w3.org/2001/XInclude"
  xmlns:exsl="http://exslt.org/common"
  extension-element-prefixes="exsl">
  
  <xsl:output method="xml"/>
  
  <xsl:template match="/">
    <map>
      <xsl:apply-templates/>
      <xsl:if test="1=2">
        <para>This is only for <emph>testing</emph> this script.</para>
      </xsl:if>
    </map>
  </xsl:template>

  <!-- don't set position=right everywhere, call a template to emit
       the attribute only for children of the root element -->

  <xsl:template match="*">
    <xsl:element name="node">
      <xsl:attribute name="POSITION">right</xsl:attribute>
      <xsl:if test="count(ancestor::*) &gt; 0">
        <xsl:attribute name="FOLDED">true</xsl:attribute>        
      </xsl:if>
      <xsl:attribute name="COLOR">#999999</xsl:attribute>
      <xsl:attribute name="TEXT">
        <xsl:value-of select="name()"/>
      </xsl:attribute>
      <xsl:apply-templates select="@*"/>
      <font NAME="Dialog" SIZE="10"/>
      <xsl:apply-templates select="*|text()"/>
    </xsl:element>
  </xsl:template>

  <xsl:template match="@*">
    <xsl:element name="node">
      <xsl:attribute name="POSITION">right</xsl:attribute>
      <xsl:attribute name="FOLDED">false</xsl:attribute>
      <xsl:attribute name="COLOR">#9999ff</xsl:attribute>
      <xsl:attribute name="TEXT">
        <xsl:value-of select="concat(name(),':')"/>
      </xsl:attribute>
      <font NAME="Dialog" SIZE="10"/>
      <xsl:element name="node">
        <xsl:attribute name="FOLDED">false</xsl:attribute>
        <xsl:attribute name="TEXT">
          <xsl:value-of select="."/>
        </xsl:attribute>
      </xsl:element>
    </xsl:element>
  </xsl:template>

  <xsl:template match="text()[normalize-space(.)]">
    <xsl:element name="node">
      <xsl:attribute name="TEXT">
        <xsl:value-of select="."/>
      </xsl:attribute>
    </xsl:element>
  </xsl:template>

</xsl:stylesheet>
