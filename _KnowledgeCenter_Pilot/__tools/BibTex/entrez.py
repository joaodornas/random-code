#!/usr/bin/env python
#
# Copyright (C) 1998-2004 Frederic Gobry
# Copyright (C) 2008-2012 W. Trevor King
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with This program.  If not, see <http://www.gnu.org/licenses/>.

"""Python interface to Entrez_ SOAP_ using the suds_ module.

Before you use this program, read the rules_.

.. _Entrez: http://www.ncbi.nlm.nih.gov/books/NBK25500/
.. _SOAP: http://www.ncbi.nlm.nih.gov/entrez/query/static/esoap_help.html
.. _suds: https://fedorahosted.org/suds/
.. _rules: http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html#UserSystemRequirements

To discover services using suds, try:

>>> print EUTILS_CLIENT  # doctest: +ELLIPSIS, +REPORT_UDIFF
<BLANKLINE>
Suds ( https://fedorahosted.org/suds/ )  version: ...  build: ...
<BLANKLINE>
Service ( eUtilsService ) tns="http://www.ncbi.nlm.nih.gov/soap/eutils/"
   Prefixes (6)
      ns0 = "http://www.ncbi.nlm.nih.gov/soap/eutils/egquery"
      ns1 = "http://www.ncbi.nlm.nih.gov/soap/eutils/einfo"
      ns2 = "http://www.ncbi.nlm.nih.gov/soap/eutils/elink"
      ns3 = "http://www.ncbi.nlm.nih.gov/soap/eutils/epost"
      ns4 = "http://www.ncbi.nlm.nih.gov/soap/eutils/esearch"
      ns5 = "http://www.ncbi.nlm.nih.gov/soap/eutils/esummary"
   Ports (1):
      (eUtilsServiceSoap)
         Methods (7):
            run_eGquery(xs:string term, xs:string tool, xs:string email, )
            run_eInfo(xs:string db, xs:string tool, xs:string email, )
            run_eLink(xs:string db, xs:string[] id, xs:string reldate, ...)
            run_ePost(xs:string db, xs:string id, xs:string WebEnv, ...)
            run_eSearch(xs:string db, xs:string term, xs:string WebEnv, ...)
            run_eSpell(xs:string db, xs:string term, xs:string tool, ...)
            run_eSummary(xs:string db, xs:string id, xs:string WebEnv, ...)
         Types (34):
            ns1:DbInfoType
            ns1:DbListType
            ...
            ns0:eGQueryResultType
<BLANKLINE>
<BLANKLINE>
"""

import logging as _logging
from xml.sax.saxutils import unescape as _unescape
import subprocess as _subprocess
import sys as _sys
import time as _time
import urllib as _urllib

import suds as _suds
from suds.client import Client as _Client
from suds.transport import TransportError as _TransportError

# Platform constants
_MSWINDOWS = _sys.platform == 'win32'
_POSIX = not _MSWINDOWS

if _POSIX:
    import os as _os
    import select as _select


__version__ = '0.2'


EUTILS_WSDL_URL = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/eutils.wsdl'
EFETCH_WSDL_URL = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/efetch_%s.wsdl'
EFETCH_PLAIN_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
NCBI_PLAIN_URL = 'http://www.ncbi.nlm.nih.gov/%s/%s'

EUTILS_CLIENT = _Client(EUTILS_WSDL_URL)

# Entrez-requested tracking information
TOOL = 'entrezpy'
EMAIL = 'wking@drexel.edu'

# Logging
LOG = _logging.getLogger(TOOL)
LOG.setLevel(_logging.WARN)
_handler = _logging.StreamHandler()
_formatter = _logging.Formatter('%(name)-8s: %(levelname)-6s %(message)s')
_handler.setFormatter(_formatter)
LOG.addHandler(_handler)
del _handler, _formatter



## Use the external bibutils package to convert to BibTeX format


class Pipe (object):
    """Simple interface for executing POSIX-style pipes.

    Based on the subprocess module.  The only complication is the
    adaptation of `subprocess.Popen._communicate` to listen to the
    stderrs of all processes involved in the pipe, as well as the
    terminal process' stdout.  There are two implementations of
    `Pipe._communicate`, one for MS Windows, and one for POSIX
    systems.  The MS Windows implementation is currently untested.

    >>> p = Pipe([['find', '/etc/'], ['grep', '^/etc/ssh$']])
    >>> p.stdout
    '/etc/ssh\\n'
    >>> p.status
    1
    >>> p.statuses
    [1, 0]
    >>> p.stderrs # doctest: +ELLIPSIS
    [...find: ...: Permission denied..., '']

    >>> p = Pipe([['cat'], ['head']], stdin='line 1\\nline 2\\nline 3\\n')
    >>> p.stdout
    'line 1\\nline 2\\nline 3\\n'
    >>> p.statuses
    [0, 0]
    >>> p.stderrs
    ['', '']
    """
    def __init__(self, cmds, stdin=None):
        if isinstance(stdin, str):
            stdin_str = stdin
            stdin = _subprocess.PIPE
        else:
            stdin_str = None

        # spawn processes
        self._procs = []
        for cmd in cmds:
            if len(self._procs) != 0:
                stdin = self._procs[-1].stdout
            LOG.debug('run command %s' % cmd)
            kwargs = {}
            if _POSIX:
                kwargs['close_fds'] = True
            try:
                self._procs.append(_subprocess.Popen(
                        cmd, stdin=stdin, stdout=_subprocess.PIPE,
                        stderr=_subprocess.PIPE, **kwargs))
            except OSError:
                LOG.error(cmd)
                raise
        self.stdout,self.stderrs = self._communicate(input=stdin_str)

        # collect process statuses
        self.statuses = []
        self.status = 0
        for proc in self._procs:
            self.statuses.append(proc.wait())
            LOG.debug('join %s (status %d)' % (proc, self.statuses[-1]))
            if self.statuses[-1] != 0:
                self.status = self.statuses[-1]

    # Code excerpted from subprocess.Popen._communicate()
    if _MSWINDOWS == True:
        def _communicate(self, input=None):
            LOG.debug('communicate with pipe')
            assert input == None, 'stdin != None not yet supported'
            # listen to each process' stderr
            threads = []
            std_X_arrays = []
            for proc in self._procs:
                stderr_array = []
                thread = Thread(target=proc._readerthread,
                                args=(proc.stderr, stderr_array))
                thread.setDaemon(True)
                thread.start()
                threads.append(thread)
                std_X_arrays.append(stderr_array)

            # also listen to the last processes stdout
            stdout_array = []
            thread = Thread(target=proc._readerthread,
                            args=(proc.stdout, stdout_array))
            thread.setDaemon(True)
            thread.start()
            threads.append(thread)
            std_X_arrays.append(stdout_array)

            # join threads as they die
            for thread in threads:
                thread.join()

            # read output from reader threads
            std_X_strings = []
            for std_X_array in std_X_arrays:
                std_X_strings.append(std_X_array[0])

            stdout = std_X_strings.pop(-1)
            stderrs = std_X_strings
            LOG.debug('pipe communication complete')
            return (stdout, stderrs)
    else:
        assert _POSIX==True, 'invalid platform'
        def _communicate(self, input=None):
            LOG.debug('communicate with pipe')
            read_set = []
            write_set = []
            read_arrays = []
            stdout = None # Return
            stderr = None # Return

            if self._procs[0].stdin:
                # Flush stdio buffer.  This might block, if the user has
                # been writing to .stdin in an uncontrolled fashion.
                self._procs[0].stdin.flush()
                if input:
                    write_set.append(self._procs[0].stdin)
                else:
                    self._procs[0].stdin.close()
            for proc in self._procs:
                read_set.append(proc.stderr)
                read_arrays.append([])
            read_set.append(self._procs[-1].stdout)
            read_arrays.append([])

            input_offset = 0
            while read_set or write_set:
                LOG.debug('select on read %s, write %s' %(read_set, write_set))
                try:
                    rlist,wlist,xlist = _select.select(read_set, write_set, [])
                except _select.error, e:
                    if e.args[0] == errno.EINTR:
                        LOG.debug('EINTR: %s' % e)
                        continue
                    raise
                LOG.debug('selected read %s, write %s, exception %s'
                          % (rlist, wlist, xlist))
                if self._procs[0].stdin in wlist:
                    # When select has indicated that the file is writable,
                    # we can write up to PIPE_BUF bytes without risk
                    # blocking.  POSIX defines PIPE_BUF >= 512
                    LOG.debug('write to stdin for process 0')
                    chunk = input[input_offset:input_offset+512]
                    bytes_written = _os.write(
                        self._procs[0].stdin.fileno(), chunk)
                    input_offset += bytes_written
                    if input_offset >= len(input):
                        self._procs[0].stdin.flush()
                        self._procs[0].stdin.close()
                        write_set.remove(self._procs[0].stdin)
                        LOG.debug('stdin complete')
                if self._procs[-1].stdout in rlist:
                    LOG.debug('read stdout for final process')
                    data = _os.read(self._procs[-1].stdout.fileno(), 1024)
                    if data == '':
                        self._procs[-1].stdout.close()
                        read_set.remove(self._procs[-1].stdout)
                        LOG.debug('stdout complete')
                    read_arrays[-1].append(data)
                for i,proc in enumerate(self._procs):
                    if proc.stderr in rlist:
                        LOG.debug('read stderr for process %i' % i)
                        data = _os.read(proc.stderr.fileno(), 1024)
                        if data == '':
                            proc.stderr.close()
                            read_set.remove(proc.stderr)
                            LOG.debug('stderr complete for process %d' % i)
                        read_arrays[i].append(data)

            # All data exchanged.  Translate lists into strings.
            read_strings = []
            for read_array in read_arrays:
                read_strings.append(''.join(read_array))

            stdout = read_strings.pop(-1)
            stderrs = read_strings
            LOG.debug('pipe communication complete')
            return (stdout, stderrs)


def medline_xml_to_bibtex(fetch_page):
    """Convert medline XML to BibTeX

    >>> xml = '\\n'.join([
    ...     '<?xml version="1.0"?>',
    ...     '<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, '
    ...     '1st January 2011//EN" "http://www.ncbi.nlm.nih.gov/entrez/query'
    ...     '/DTD/pubmed_110101.dtd">',
    ...     '<PubmedArticleSet>',
    ...     ' <PubmedArticle>',
    ...     '  <MedlineCitation Owner="NLM" Status="MEDLINE">',
    ...     '   <PMID Version="1">20004685</PMID>',
    ...     '   <Article PubModel="Print-Electronic">',
    ...     '    <Journal>',
    ...     '     <ISSN IssnType="Electronic">1879-0003</ISSN>',
    ...     '     <JournalIssue CitedMedium="Internet">',
    ...     '      <Volume>46</Volume><Issue>2</Issue>',
    ...     '      <PubDate>',
    ...     '       <Year>2010</Year><Month>Mar</Month><Day>1</Day>',
    ...     '      </PubDate>',
    ...     '     </JournalIssue>',
    ...     '    </Journal>',
    ...     '    <ArticleTitle>Monte Carlo simulation of mechanical unfolding '
    ...          'of proteins based on a simple two-state model.'
    ...          '</ArticleTitle>',
    ...     '    <Pagination><MedlinePgn>159-66</MedlinePgn></Pagination>',
    ...     '    <AuthorList CompleteYN="Y">',
    ...     '     <Author ValidYN="Y">',
    ...     '      <LastName>King</LastName>',
    ...     '      <ForeName>William T</ForeName>',
    ...     '      <Initials>WT</Initials>',
    ...     '     </Author>',
    ...     '     <Author ValidYN="Y">',
    ...     '      <LastName>Su</LastName>',
    ...     '      <ForeName>Meihong</ForeName>',
    ...     '      <Initials>M</Initials>',
    ...     '     </Author>',
    ...     '     <Author ValidYN="Y">',
    ...     '      <LastName>Yang</LastName>',
    ...     '      <ForeName>Guoliang</ForeName>',
    ...     '      <Initials>G</Initials>',
    ...     '     </Author>',
    ...     '    </AuthorList>',
    ...     '    <MedlineJournalInfo>',
    ...     '     <MedlineTA>Int J Biol Macromol</MedlineTA>',
    ...     '    </MedlineJournalInfo>',
    ...     '   </Article>',
    ...     '   <MedlineJournalInfo>',
    ...     '    <MedlineTA>Int J Biol Macromol</MedlineTA>',
    ...     '   </MedlineJournalInfo>',
    ...     '  </MedlineCitation>',
    ...     '  <PubmedData>',
    ...     '   <ArticleIdList>',
    ...     '    <ArticleId IdType="doi">10.1016/j.ijbiomac.2009.12.001'
    ...          '</ArticleId>',
    ...     '   </ArticleIdList>',
    ...     '  </PubmedData>',
    ...     ' </PubmedArticle>',
    ...     '</PubmedArticleSet>',
    ...     ])
    >>> print medline_xml_to_bibtex(xml)  # doctest: +REPORT_UDIFF
    @Article{King2010,
      author =       "William T. King and Meihong Su and Guoliang Yang",
      title =        "Monte Carlo simulation of mechanical unfolding of
                     proteins based on a simple two-state model.",
      journal =      "Int J Biol Macromol",
      year =         "2010",
      month =        mar,
      day =          "01",
      volume =       "46",
      number =       "2",
      pages =        "159--166",
      ISSN =         "1879-0003",
      doi =          "10.1016/j.ijbiomac.2009.12.001",
      URL =          "http://www.ncbi.nlm.nih.gov/pubmed/20004685",
    }
    <BLANKLINE>
    """
    LOG.info('convert medline XML to BibTeX')
    LOG.debug('convert from\n%s' % fetch_page)
    p = Pipe(cmds=[['med2xml'], ['xml2bib', '-fc'], ['bibclean']],
             stdin=fetch_page)
    LOG.debug('converted to\n%s' % p.stdout)
    return p.stdout


if __name__ == '__main__':
    from optparse import OptionParser

    usage_string = '\n'.join([
            '',
            '  %prog [options] SEARCH_TERM'
             '       (print medline xml matching search)',
            '| %prog -l [options] SEARCH_TERM'
             '    (print links to entries matching search)',
            '| %prog -L [-f FILE]                (list databases)',
            '| %prog -X [-d DATABASE] [-F FIELD] [-f FILE]'
             '  (list fields in a database, or details on a single field)',
            '',
            '2008-2011, W. Trevor King.',
            '',
            'See the docstrings in %prog or',
            ' http://www.ncbi.nlm.nih.gov/books/NBK3837/',
            ' http://www.ncbi.nlm.nih.gov/entrez/query/static/'
             'eutils_help.html',
            ' http://www.ncbi.nlm.nih.gov/entrez/query/static/'
             'eutils_help.html#UserSystemRequirements',
            ' http://www.ncbi.nlm.nih.gov/corehtml/query/static/'
             'einfo_help.html',
            ' http://www.ncbi.nlm.nih.gov/corehtml/query/static/'
             'esearch_help.html',
            ' http://www.ncbi.nlm.nih.gov/corehtml/query/static/'
             'efetch_help.html',
            ' http://www.ncbi.nlm.nih.gov/corehtml/query/static/'
             'elink_help.html',
            'for more details.'
            ])

    parser = OptionParser(
        usage=usage_string, version='%%prog %s' % __version__)

    # Explaination by Jerry Stratton, http://www.hoboes.com/Mimsy/?ART=511
    # "
    # metavar is the name used in the help for that options required
    # text, and dest is the name of the property you'll use to access
    # the value of that option.
    # "

    parser.add_option('-d', '--database', dest='database',
                      help="Search DATABASE (default '%default')",
                      type='string', metavar='DATABASE', default='pubmed')
    parser.add_option('-f', '--file', dest='filename',
                      help='write output to FILE (default stdout)',
                      type='string', metavar='FILE')
    parser.add_option('-v', '--verbose', dest='verbose', action='count',
                      help=('Print minimal debugging information.  Use twice '
                            'to get lots of debugging info.'),
                      default=0)

    # mode control options
    mode = 'search'
    def set_mode(option, opt_str, value, parser):
        global mode
        long_option = option.get_opt_string()
        if long_option == '--list-mode':
            mode = 'list'
        elif long_option == '--explain-mode':
            mode = 'explain'

    parser.add_option('-L', '--list-mode', callback=set_mode,
                      help='Run in list mode', action='callback')
    parser.add_option('-X', '--explain-mode', callback=set_mode,
                      help='Run in explain mode', action='callback')

    # search-fetch-xml-to-? options
    output = 'bibtex'
    def set_output(option, opt_str, value, parser):
        global output
        long_option = option.get_opt_string()
        if long_option == '--output-xml':
            output = 'medline'
        if long_option == '--output-bibtex':
            output = 'bibtex'
        if long_option == '--output-link':
            output = 'link'
    parser.add_option('-x', '--output-xml', callback=set_output,
                      help='Output search results as Medline XML',
                      action='callback')
    parser.add_option('-b', '--output-bibtex', callback=set_output,
                      help='Output search results as BibTeX',
                      action='callback')
    parser.add_option('-F', '--field', dest='field',
                      help='Limit SEARCH_TERM to FIELD',
                      type='string', metavar='FIELD')
    parser.add_option('-r', '--reldate', dest='reldate',
                      help='Limit search to dates within DAYS of today',
                      type='string', metavar='DAYS')
    parser.add_option('--mindate', dest='mindate',
                      help=('Limit search to date after MINDATE '
                            "(e.g. '2001/1/1' or '2002')"),
                      type='string', metavar='MINDATE')
    parser.add_option('--maxdate', dest='maxdate',
                      help=('Limit search to date after MAXDATE '
                            "(e.g. '2001/1/1' or '2002')"),
                      type='string', metavar='MAXDATE')
    parser.add_option('-t', '--datetype', dest='datetype',
                      help=("Select field to apply date limits to "
                            "(e.g. 'edat' for Entrez date)"),
                      type='string', metavar='DATETYPE')
    parser.add_option('-m', '--retmax', dest='retmax',
                      help=('Return at most RETMAX items from a successful '
                            'search (default %default)'),
                      type='int', metavar='RETMAX', default=20)
    parser.add_option('-s', '--retstart', dest='retstart',
                      help=('Index of first returned search item from a '
                            'successful search (default %default)'),
                      type='int', metavar='RETSTART', default=0)
    parser.add_option('-V', '--validate', dest='validate', action='store_true',
                      help=('Check that FIELD and field tags in SEARCH_TERM '
                            'are valid for DB'),
                      default=False)

    # output link options
    parser.add_option('-l', '--output-link', callback=set_output,
                      help='Output a link (instead of xml citations).',
                      action='callback')
    parser.add_option('-c', '--link-cmd', dest='link_cmd',
                      help='Select link output',
                      type='string', metavar='LINK_CMD')
    parser.add_option('-T', '--link-term', dest='link_term',
                      help='Limit links to those matching LINK_TERM',
                      type='string', metavar='LINK_TERM')
    parser.add_option('-D', '--from-database', dest='dbfrom',
                      help='Limit links to those from FROMDATABASE)',
                      type='string', metavar='FROMDATABASE')
    parser.add_option('-n', '--link-name', dest='linkname',
                      help='Limit links to a specific neighbor',
                      type='string', metavar='LINKNAME')

    (options, args) = parser.parse_args()
    parser.destroy()

    # open the output file if specified
    if options.filename == None:
        outfile = _sys.stdout
    else:
        outfile = file(options.filename, 'w')

    if options.verbose == 1:
        LOG.setLevel(_logging.INFO)
    elif options.verbose > 1:
        LOG.setLevel(_logging.DEBUG)

    LOG.debug('operating in %s mode' % mode)

    if mode == 'list':
        outfile.write('# available databases:\n')
        LOG.info('run eInfo to get list of databases')
        q = EUTILS_CLIENT.service.run_eInfo(tool=TOOL, email=EMAIL)
        if hasattr(q, 'ERROR'):
            raise Exception(q.ERROR)

        for db in q.DbList.DbName:
            outfile.write('%s\n' % db)

    elif mode == 'explain':
        LOG.info('run eInfo on %s' % options.database)
        q = EUTILS_CLIENT.service.run_eInfo(
            db=options.database, tool=TOOL, email=EMAIL)
        if hasattr(q, 'ERROR'):
            raise Exception(q.ERROR)

        if options.field:  # print specific info about this field
            outfile.write(
                'field %s in %s:\n' % (options.field, options.database))
            fields = dict(
                [(field.Name, field) for field in q.DbInfo.FieldList.Field])
            field = fields[options.field]
            attributes = sorted(
                [(a, getattr(field, a)) for a in dir(field)
                 if not a.startswith('_')])
            field_size = [0]
            for attribute,value in attributes:
                if len(attribute) > field_size[0]:
                    field_size[0] = len(attribute)
            for attribute,value in attributes:
                outfile.write(
                    '%*.*s\t%s\n'
                    % (field_size[0], field_size[0], attribute, value))
        else:  # print general info
            outfile.write('database: %s\n' % q.DbInfo.DbName)
            outfile.write('description: %s\n' % q.DbInfo.Description)
            outfile.write('available fields:\n')
            field_size = [0,0]
            for field in q.DbInfo.FieldList.Field:
                if len(field.Name) > field_size[0]:
                    field_size[0] = len(field.Name)
                if len(field.FullName) > field_size[1]:
                    field_size[1] = len(field.FullName)
            for field in q.DbInfo.FieldList.Field:
                outfile.write(
                    '%*.*s\t%-*.*s\t%s\n'
                    % (field_size[0], field_size[0], field.Name,
                       field_size[1], field_size[1], field.FullName,
                       field.Description))

    elif mode == 'search':
        search_term = args[0]
        LOG.debug('output %s' % output)

        if options.mindate and not options.maxdate:
            options.maxdate = _time.strftime('%Y/%M/%d')
            LOG.info('fill in maximum date: %s' % options.maxdate)
        elif options.maxdate and not options.mindate:
            options.mindate = '0'
            LOG.info('fill in minimum date: %s' % options.mindate)

        LOG.info('run eEsearch on %s' % options.database)
        q = EUTILS_CLIENT.service.run_eSearch(
            db=options.database, term=search_term, tool=TOOL, email=EMAIL,
            field=options.field, reldate=options.reldate,
            mindate=options.mindate, maxdate=options.maxdate,
            datetype=options.datetype, 
            RetStart=options.retstart, RetMax=options.retmax,
            #sort=)
            )
        if hasattr(q, 'ERROR'):
            raise Exception(q.ERROR)
        if hasattr(q.IdList, 'Id'):
            ret = int(len(q.IdList.Id))
        else:
            ret = 0
        LOG.info('search returned %d of %d items' % (ret, int(q.Count)))

        if ret > 0:
            if output in ['medline', 'bibtex']:
                e = None
                try:
                    efetch_client = _Client(EFETCH_WSDL_URL % options.database)
                except _TransportError, e:
                    if e.httpcode != 404:
                        raise
                    LOG.warn(str(e))
                if e:  # Fallback to straight URL fetch
                    params = {
                        'id': ','.join(q.IdList.Id),
                        'tool': TOOL,
                        'email': EMAIL,
                        'db': options.database,
                        'report': 'xml',
                        }
                    url = '%s?%s' % (
                        EFETCH_PLAIN_URL, _urllib.urlencode(params))
                    LOG.info('fallback to non-SOAP eFetch request: %s' % url)
                    f = _urllib.urlopen(url)
                    xml = f.read()
                    f.close()
                    # Remove wrapping HTML and unescape XML
                    #LOG.debug('raw data:\n%s' % xml)
                    xml = xml.split('<pre>', 1)[-1]
                    xml = xml.split('</pre>', 1)[0]
                    xml = _unescape(xml, {'&quot;': '"'})
                    #LOG.debug('xml data:\n%s' % xml)
                    if not xml.strip(): # 
                        urls = [NCBI_PLAIN_URL % (options.database, id)
                                for id in q.IdList.Id]
                        LOG.warn(
                            'no meaningful output; try:\n%s' % '\n'.join(urls))
                else:  # Use SOAP eFetch
                    LOG.info('run eFetch on %s' % options.database)
                    f = efetch_client.service.run_eFetch(
                        id=','.join(q.IdList.Id), tool=TOOL, email=EMAIL)
                    if hasattr(f, 'ERROR'):
                        raise Exception(f.ERROR)
                    xml = efetch_client.last_received()

            if output is None:
                pass  # we're bailing
            elif output == 'medline':
                outfile.write(str(xml).rstrip()+'\n')
            elif output == 'bibtex':
                outfile.write(medline_xml_to_bibtex(str(xml)))
            elif output == 'link':
                LOG.info('run eLink on %s' % options.database)
                f = EUTILS_CLIENT.service.run_eLink(
                    db=options.database, id=','.join(q.IdList.Id),
                    #reldate=, mindate=, maxdate=, datetype=,
                    term=options.link_term, dbfrom=options.dbfrom,
                    linkname=options.linkname, cmd=options.link_cmd,
                    tool=TOOL, email=EMAIL)
                outfile.write(str(EUTILS_CLIENT.last_received()).rstrip()+'\n')
            else:
                raise KeyError(output)

    if options.filename != None:
        outfile.close()
