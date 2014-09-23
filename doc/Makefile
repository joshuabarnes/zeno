# ~/zeno/doc/Makefile: generate documentation for zeno programs.

build_all:
	$(MAKE) -s start_table TITLE="Zeno Programs"
	$(MAKE) -s build_header HEADER="N-Body Simulation"
	$(MAKE) -s build_table build_page NAME=treecode
	$(MAKE) -s build_table build_page NAME=treecode_x
	$(MAKE) -s build_table build_page NAME=treepeek
	$(MAKE) -s build_header HEADER="SPH Simulation"
	$(MAKE) -s build_table build_page NAME=sphcode_e
	$(MAKE) -s build_table build_page NAME=sphcode_ea
	$(MAKE) -s build_table build_page NAME=sphcode_u
	$(MAKE) -s build_table build_page NAME=sphcode_ui
	$(MAKE) -s build_table build_page NAME=sphcode_ux
	$(MAKE) -s build_table build_page NAME=sphcode_ur
	$(MAKE) -s build_table build_page NAME=sphcode_uo
	$(MAKE) -s build_table build_page NAME=sphcode_ud
	$(MAKE) -s build_table build_page NAME=sphcode_uc
	$(MAKE) -s build_table build_page NAME=sphcode_us
	$(MAKE) -s build_table build_page NAME=sphcode_uis
	$(MAKE) -s build_table build_page NAME=sphcode_uir
	$(MAKE) -s build_header HEADER="Snapshot Construction"
	$(MAKE) -s build_table build_page NAME=kepsnap
	$(MAKE) -s build_table build_page NAME=makesnap
	$(MAKE) -s build_table build_page NAME=plumsnap
	$(MAKE) -s build_table build_page NAME=polysnap
	$(MAKE) -s build_table build_page NAME=snaptestdisk
	$(MAKE) -s build_table build_page NAME=snapcons
	$(MAKE) -s build_table build_page NAME=snapstack
	$(MAKE) -s build_header HEADER="General Spherical Profile Generation"
	$(MAKE) -s build_table build_page NAME=expdgsp
	$(MAKE) -s build_table build_page NAME=gammagsp
	$(MAKE) -s build_table build_page NAME=halogsp
	$(MAKE) -s build_table build_page NAME=isothgsp
	$(MAKE) -s build_table build_page NAME=plumgsp
	$(MAKE) -s build_table build_page NAME=polygsp
	$(MAKE) -s build_table build_page NAME=snapgsp
	$(MAKE) -s build_table build_page NAME=tablegsp
	$(MAKE) -s build_header HEADER="General Spherical Profile Utilities"
	$(MAKE) -s build_table build_page NAME=gspadd
	$(MAKE) -s build_table build_page NAME=gspeval
	$(MAKE) -s build_table build_page NAME=gsplist
	$(MAKE) -s build_table build_page NAME=gspsmooth
	$(MAKE) -s build_table build_page NAME=gspsmoothfit
	$(MAKE) -s build_table build_page NAME=gsptrun
	$(MAKE) -s build_header HEADER="General Spherical Profile Realization"
	$(MAKE) -s build_table build_page NAME=gspdisk
	$(MAKE) -s build_table build_page NAME=gspmodel
	$(MAKE) -s build_table build_page NAME=gsprealize
	$(MAKE) -s build_table build_page NAME=gspsphere
	$(MAKE) -s build_table build_page NAME=gsptestdisk
	$(MAKE) -s build_header HEADER="Snapshot Assignment"	
	$(MAKE) -s build_table build_page NAME=snapaux
	$(MAKE) -s build_table build_page NAME=snapkey
	$(MAKE) -s build_table build_page NAME=snapset
	$(MAKE) -s build_table build_page NAME=buildmap
	$(MAKE) -s build_header HEADER="Snapshot Transformation"	
	$(MAKE) -s build_table build_page NAME=snapavg
	$(MAKE) -s build_table build_page NAME=snapcenter
	$(MAKE) -s build_table build_page NAME=snapcollect
	$(MAKE) -s build_table build_page NAME=snapforce
	$(MAKE) -s build_table build_page NAME=snapobject
	$(MAKE) -s build_table build_page NAME=snaprect
	$(MAKE) -s build_table build_page NAME=snaprotate
	$(MAKE) -s build_table build_page NAME=snapsift
	$(MAKE) -s build_table build_page NAME=snapsort
	$(MAKE) -s build_table build_page NAME=snapsplit
	$(MAKE) -s build_table build_page NAME=snapstream
	$(MAKE) -s build_table build_page NAME=snaptrak
	$(MAKE) -s build_table build_page NAME=snapturn
	$(MAKE) -s build_header HEADER="Snapshot Analysis"	
	$(MAKE) -s build_table build_page NAME=snapcm
	$(MAKE) -s build_table build_page NAME=snaphist
	$(MAKE) -s build_table build_page NAME=snapshape
	$(MAKE) -s build_table build_page NAME=snapstat
	$(MAKE) -s build_table build_page NAME=snaprholist
	$(MAKE) -s build_header HEADER="Snapshot Output"	
	$(MAKE) -s build_table build_page NAME=snapascii
	$(MAKE) -s build_table build_page NAME=snaplist
	$(MAKE) -s build_table build_page NAME=snaptipsy
	$(MAKE) -s build_header HEADER="Snapshot Visuialization"	
	$(MAKE) -s build_table build_page NAME=snapppm
	$(MAKE) -s build_table build_page NAME=snapsmooth
	$(MAKE) -s build_table build_page NAME=snapview
	$(MAKE) -s build_table build_page NAME=viewfile
	$(MAKE) -s build_header HEADER="Structured File Utilities"
	$(MAKE) -s build_table build_page NAME=csf
	$(MAKE) -s build_table build_page NAME=dsf
	$(MAKE) -s build_table build_page NAME=gsf
	$(MAKE) -s build_table build_page NAME=rsf
	$(MAKE) -s build_table build_page NAME=ssf
	$(MAKE) -s build_table build_page NAME=tsf
	$(MAKE) -s end_table

build_page:
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\""                 >  $(NAME).html
	echo "\"http://www.w3.org/TR/html4/loose.dtd\">"                                        >> $(NAME).html
	echo "<html>"                                                                           >> $(NAME).html
	echo "<head>"                                                                           >> $(NAME).html
	echo "  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">"   >> $(NAME).html
	echo "  <title>$(NAME)</title>"                                                         >> $(NAME).html
	echo "  <style type=\"text/css\">"                                                      >> $(NAME).html
	echo "    a {text-decoration: none}"                                                    >> $(NAME).html
	echo "  </style>"                                                                       >> $(NAME).html
	echo "</head>"                                                                          >> $(NAME).html
	echo "<body bgcolor=\"white\">"                                                         >> $(NAME).html
	echo "<h1>$(NAME)</h1>"                                                                 >> $(NAME).html
	echo "<hr>"                                                                             >> $(NAME).html
	echo "<pre>"                                                                            >> $(NAME).html
	$(NAME) -explain                                                                        >> $(NAME).html
	echo "</pre>"                                                                           >> $(NAME).html
	echo "<hr>"                                                                             >> $(NAME).html
	echo "<table width=\"100%\">"                                                           >> $(NAME).html
	echo "  <tr valign=\"top\">"                                                            >> $(NAME).html
	echo "    <td align=\"left\">"                                                          >> $(NAME).html
	echo "      <address>"                                                                  >> $(NAME).html
	echo "        <a href=\"http://www.ifa.hawaii.edu/~barnes/\">Joshua E. Barnes</a>"      >> $(NAME).html
	echo "        &nbsp;&nbsp;&nbsp;&nbsp;"                                                 >> $(NAME).html
	echo "        (<tt>barnes at ifa.hawaii.edu</tt>)"                                      >> $(NAME).html
	echo "      </address>"                                                                 >> $(NAME).html
	echo "      Created: "                                                                  >> $(NAME).html
	date +"      %e %B %Y"                                                                  >> $(NAME).html
	echo "      <br>"                                                                       >> $(NAME).html
	echo "      <tt>http://www.ifa.hawaii.edu/~barnes/zeno/$(NAME).html</tt>"               >> $(NAME).html
	echo "    </td>"                                                                        >> $(NAME).html
	echo "    <td align=\"right\">"                                                         >> $(NAME).html
	echo "      <a href=\"http://validator.w3.org/check?uri=referer\"><img border=\"0\""    >> $(NAME).html
	echo "        src=\"http://www.w3.org/Icons/valid-html401\""                            >> $(NAME).html
	echo "        alt=\"Valid HTML 4.01!\" height=\"31\" width=\"88\"></a>"                 >> $(NAME).html
	echo "    </td>"                                                                        >> $(NAME).html
	echo "  </tr>"                                                                          >> $(NAME).html
	echo "</table>"                                                                         >> $(NAME).html
	echo "</body>"                                                                          >> $(NAME).html
	echo "</html>"                                                                          >> $(NAME).html

TABLE = Programs.html

start_table: ProgramText.txt
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\""                 >  $(TABLE)
	echo "\"http://www.w3.org/TR/html4/loose.dtd\">"                                        >> $(TABLE)
	echo "<html>"                                                                           >> $(TABLE)
	echo "<head>"                                                                           >> $(TABLE)
	echo "  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">"   >> $(TABLE)
	echo "  <title>$(TITLE)</title>"                                                        >> $(TABLE)
	echo "  <style type=\"text/css\">"                                                      >> $(TABLE)
	echo "    a {text-decoration: none}"                                                    >> $(TABLE)
	echo "  </style>"                                                                       >> $(TABLE)
	echo "</head>"                                                                          >> $(TABLE)
	echo "<body bgcolor=\"white\">"                                                         >> $(TABLE)
	echo "<h1>$(TITLE)</h1>"                                                                >> $(TABLE)
	echo "<hr>"                                                                             >> $(TABLE)
	cat ProgramText.txt                                                                     >> $(TABLE)
	echo "<table border=\"0\" cellpadding=\"0\">"                                           >> $(TABLE)

build_header:
	echo "<tr> <td colspan=\"4\"><b>$(HEADER)</b></td> </tr>"                               >> $(TABLE)

build_table:
	$(NAME) -help | \
	  head -1 | \
	    gawk '{ \
	            printf("<tr> <td width=\"5%\">&nbsp;</td> <td><a href=\"%s.html\">%s</a></td>" \
	                   "<td width=\"5%\">&nbsp;</td> <td>%s</td> </tr>\n", \
	                   $$1, $$1, substr($$0, 35)); \
	          }'                                                                            >> $(TABLE)

end_table:
	echo "</table>"                                                                         >> $(TABLE)
	echo "<hr>"                                                                             >> $(TABLE)
	echo "<table width=\"100%\">"                                                           >> $(TABLE)
	echo "  <tr valign=\"top\">"                                                            >> $(TABLE)
	echo "    <td align=\"left\">"                                                          >> $(TABLE)
	echo "      <address>"                                                                  >> $(TABLE)
	echo "        <a href=\"http://www.ifa.hawaii.edu/~barnes/\">Joshua E. Barnes</a>"      >> $(TABLE)
	echo "        &nbsp;&nbsp;&nbsp;&nbsp;"                                                 >> $(TABLE)
	echo "        (<tt>barnes at ifa.hawaii.edu</tt>)"                                      >> $(TABLE)
	echo "      </address>"                                                                 >> $(TABLE)
	echo "      Created: "                                                                  >> $(TABLE)
	date +"      %e %B %Y"                                                                  >> $(TABLE)
	echo "      <br>"                                                                       >> $(TABLE)
	echo "      <tt>http://www.ifa.hawaii.edu/~barnes/zeno/$(TABLE)</tt>"                   >> $(TABLE)
	echo "    </td>"                                                                        >> $(TABLE)
	echo "    <td align=\"right\">"                                                         >> $(TABLE)
	echo "      <a href=\"http://validator.w3.org/check?uri=referer\"><img border=\"0\""    >> $(TABLE)
	echo "        src=\"http://www.w3.org/Icons/valid-html401\""                            >> $(TABLE)
	echo "        alt=\"Valid HTML 4.01!\" height=\"31\" width=\"88\"></a>"                 >> $(TABLE)
	echo "    </td>"                                                                        >> $(TABLE)
	echo "  </tr>"                                                                          >> $(TABLE)
	echo "</table>"                                                                         >> $(TABLE)
	echo "</body>"                                                                          >> $(TABLE)
	echo "</html>"                                                                          >> $(TABLE)
