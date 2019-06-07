<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>asg.hh</name>
    <path>/Users/symes/Applications/trip/rvl/tsopt/include/</path>
    <filename>asg_8hh</filename>
    <includes id="GridSpace_8hh" name="GridSpace.hh" local="yes" imported="no">GridSpace.hh</includes>
    <class kind="class">RVL::ASGaux</class>
    <class kind="class">RVL::ASGapplyFO</class>
    <class kind="class">RVL::ASGapplyAdjFO</class>
    <class kind="class">RVL::ASGapplyTangentFO</class>
    <class kind="class">RVL::ASGapplyAdjTangentFO</class>
    <class kind="class">RVL::ASGStep</class>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>float *</type>
      <name>sgcoeffs</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a1f19b1c5eaccb7b96b056e291d7cbdd4</anchor>
      <arglist>(int k)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Grid * &gt;</type>
      <name>make_asg_ctrllist</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a9587e9bc6f5a0a3d5cc90d487741ad9d</anchor>
      <arglist>(Grid const &amp;phys, int maxoff)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Grid * &gt;</type>
      <name>make_asg_statelist</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a78994ac742904536267943832863507a</anchor>
      <arglist>(Grid const &amp;phys, int maxoff)</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; std::string, int &gt;</type>
      <name>make_asg_indices</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a6a4ebe8273d42505edd4481f1974478b</anchor>
      <arglist>(int dim)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_pmlaxis</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a6b6b057b9b6bef3bf0158af1238ce840</anchor>
      <arglist>(int n0, int nl, int nr, float amp, float dt, int gtype, float **ep, float **epp)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>asgfcns.hh</name>
    <path>/Users/symes/Applications/trip/rvl/tsopt/include/</path>
    <filename>asgfcns_8hh</filename>
    <member kind="function">
      <type>void</type>
      <name>asg_pstep2d</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a6a04b24116b1ca73b8d3117f773cc70e</anchor>
      <arglist>(float **restrict bulk, float **restrict p0, float **restrict p1, float **restrict v0, float **restrict v1, float **restrict ep, float **restrict epp, float *restrict sdiv, const int *gsc, const int *gec, const int *lbc, const int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_pstep2d_d</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a1d1c0db36a1dc5d59edbcff1bf3d15f5</anchor>
      <arglist>(float **restrict bulk, float **restrict dbulk, float **restrict p0, float **restrict dp0, float **restrict p1, float **restrict dp1, float **restrict v0, float **restrict dv0, float **restrict v1, float **restrict dv1, float **restrict ep, float **restrict epp, float *restrict sdiv, float *restrict dsdiv, int *gsc, int *gec, int *lbc, int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_pstep2d_b</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a14007148ea16e42136aee655e0bbfc48</anchor>
      <arglist>(float **restrict bulk, float **restrict dbulk, float **restrict p0, float **restrict dp0, float **restrict p1, float **restrict dp1, float **restrict v0, float **restrict dv0, float **restrict v1, float **restrict dv1, float **restrict ep, float **restrict epp, float *restrict sdiv, float *restrict dsdiv, int *gsc, int *gec, int *lbc, int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_pstep3d</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a5d889981ec28b390cd2f4d86ef7a5270</anchor>
      <arglist>(float ***restrict bulk, float ***restrict p0, float ***restrict p1, float ***restrict p2, float ***restrict v0, float ***restrict v1, float ***restrict v2, float **restrict ep, float **restrict epp, float *restrict sdiv, const int *gsc, const int *gec, const int *lbc, const int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_vstep2d</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a4c92b6c27ca740214e7cb2c091f29ca2</anchor>
      <arglist>(float **restrict buoy, float **restrict p0, float **restrict p1, float **restrict v0, float **restrict v1, float **restrict ev, float **restrict evp, float **restrict gradp, const int *gsc_v0, const int *gec_v0, const int *gsc_v1, const int *gec_v1, const int *lbc, const int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_vstep2d_d</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a826938a4a381ade1ba4336757728318f</anchor>
      <arglist>(float **restrict buoy, float **restrict dbuoy, float **restrict p0, float **restrict dp0, float **restrict p1, float **restrict dp1, float **restrict v0, float **restrict dv0, float **restrict v1, float **restrict dv1, float **restrict ev, float **restrict evp, float **restrict gradp, float **restrict dgradp, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1, int *lbc, int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_vstep2d_b</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a44bf6a6c7e4e8ec3e6860f3a29b4d9a7</anchor>
      <arglist>(float **restrict buoy, float **restrict dbuoy, float **restrict p0, float **restrict dp0, float **restrict p1, float **restrict dp1, float **restrict v0, float **restrict dv0, float **restrict v1, float **restrict dv1, float **restrict ep, float **restrict epp, float **restrict gradp, float **restrict dgradp, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1, int *lbc, int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_vstep3d</name>
      <anchorfile>asgfcns_8hh.html</anchorfile>
      <anchor>a5f51eb517497a4c0feed48ae7990fdf9</anchor>
      <arglist>(float ***restrict buoy, float ***restrict p0, float ***restrict p1, float ***restrict p2, float ***restrict v0, float ***restrict v1, float ***restrict v2, float **restrict ev, float **restrict evp, float **restrict gradp, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1, int *gsc_v2, int *gec_v2, int *lbc, int *rbc, int maxoff, float **restrict c)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/symes/Applications/trip/rvl/tsopt/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>GridSample.hh</name>
    <path>/Users/symes/Applications/trip/rvl/tsopt/include/</path>
    <filename>GridSample_8hh</filename>
    <includes id="GridSpace_8hh" name="GridSpace.hh" local="yes" imported="no">GridSpace.hh</includes>
    <class kind="class">RVL::GridCopyOverlapFO</class>
    <class kind="class">RVL::GridExtendFO</class>
    <class kind="class">RVL::GridtoTSOp</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>GridSpace.hh</name>
    <path>/Users/symes/Applications/trip/rvl/tsopt/include/</path>
    <filename>GridSpace_8hh</filename>
    <includes id="TSOp_8hh" name="TSOp.hh" local="yes" imported="no">TSOp.hh</includes>
    <class kind="class">RVL::Axis</class>
    <class kind="class">RVL::Grid</class>
    <class kind="class">RVL::GridDC</class>
    <class kind="class">RVL::GridDCFactory</class>
    <class kind="class">RVL::GridSpace</class>
    <class kind="class">RVL::GridDomain</class>
    <class kind="class">RVL::GridDCIOFO</class>
    <namespace>RVL</namespace>
    <member kind="define">
      <type>#define</type>
      <name>TOLFAC</name>
      <anchorfile>GridSpace_8hh.html</anchorfile>
      <anchor>a6f1d333e0290f41444b8c69549106723</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>GRID_MAXDIM</name>
      <anchorfile>GridSpace_8hh.html</anchorfile>
      <anchor>a967d3b793af4ea281227245dba46596e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>EXTINT</name>
      <anchorfile>GridSpace_8hh.html</anchorfile>
      <anchor>a8cc44522d2f5b706b3af5b839c0ba3ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>areCompatible</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a253d4433eb023895798e937d483a50f8</anchor>
      <arglist>(Axis const &amp;a, Axis const &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>areCompatible</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ad1d907ab39aa73cf5c50aa7f450b7ad2</anchor>
      <arglist>(Grid const &amp;g, Grid const &amp;h)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>writeMeta&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a153d53575b3e91a4d34c9cf4a4782f36</anchor>
      <arglist>(Grid const &amp;g, ostream &amp;e)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataSize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa35397ec153fa0b4c587e9186c0a96f7</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMetaSize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a8f67c9e5d41d97dbaa6aeccde5c8ca28</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>serialize&lt; Axis &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aec2738851d7d7dded515e4a99db99cbb</anchor>
      <arglist>(Axis const &amp;a, size_t &amp;len)</arglist>
    </member>
    <member kind="function">
      <type>Axis *</type>
      <name>deserialize&lt; Axis &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa278468793367e9344c5e43af61f1f1c</anchor>
      <arglist>(char *cbuf, size_t len)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>serialize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ae07fa7da4d6b2bbbb58b6705d9d43c0b</anchor>
      <arglist>(Grid const &amp;g, size_t &amp;len)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>deserialize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ae0bf0ea161cae5d14063865df8a0a357</anchor>
      <arglist>(char *cbuf, size_t len)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Grid &gt;</type>
      <name>make_padded_grid</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a434e8bfa1dc2bb99e22acdf49ab3265f</anchor>
      <arglist>(Grid const &amp;phys, std::vector&lt; int &gt; nlsloc, std::vector&lt; int &gt; nrsloc)</arglist>
    </member>
    <member kind="function">
      <type>T **</type>
      <name>dimView</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aef62f8bbac064b68be274998dc02857d</anchor>
      <arglist>(T *x, size_t n, size_t chunklen)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Grid &gt;</type>
      <name>GridfromFile</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a167096aab0f1e00c5fdb9272ed2f5727</anchor>
      <arglist>(std::string fname)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>TSOp.hh</name>
    <path>/Users/symes/Applications/trip/rvl/tsopt/include/</path>
    <filename>TSOp_8hh</filename>
    <class kind="class">RVL::TSSample</class>
    <class kind="class">RVL::TSFO</class>
    <class kind="class">RVL::TSDC</class>
    <class kind="class">RVL::TSSpace</class>
    <class kind="class">RVL::TimeStepOp</class>
    <class kind="class">RVL::LinRestrictTSStep</class>
    <class kind="class">RVL::TanRestrictTSStep</class>
    <class kind="class">RVL::TSStep</class>
    <class kind="class">RVL::LinRestrictTSStep</class>
    <class kind="class">RVL::TanRestrictTSStep</class>
    <class kind="class">RVL::TimeStepOp</class>
    <class kind="class">RVL::TimeStepOpAllCache</class>
    <namespace>RVL</namespace>
    <member kind="define">
      <type>#define</type>
      <name>TOL</name>
      <anchorfile>TSOp_8hh.html</anchorfile>
      <anchor>a156b862ebf6d213f5da19b9e3ccb779e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASGapplyAdjFO</name>
    <filename>classRVL_1_1ASGapplyAdjFO.html</filename>
    <base>RVL::TSFO</base>
    <member kind="function">
      <type></type>
      <name>ASGapplyAdjFO</name>
      <anchorfile>classRVL_1_1ASGapplyAdjFO.html</anchorfile>
      <anchor>ab12d2b31c18b7d64e0448f2a4c169cd3</anchor>
      <arglist>(ASGaux const &amp;_aux, bool const &amp;_integerstep)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ASGapplyAdjFO</name>
      <anchorfile>classRVL_1_1ASGapplyAdjFO.html</anchorfile>
      <anchor>aa14e6fbee841c8df881b4a62002a9780</anchor>
      <arglist>(ASGapplyAdjFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASGapplyAdjFO</name>
      <anchorfile>classRVL_1_1ASGapplyAdjFO.html</anchorfile>
      <anchor>a2ba8536dd299e73c1c6b1f318a5005db</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ASGapplyAdjFO.html</anchorfile>
      <anchor>a377e8256b9d56ec03abb89909d8dcd1c</anchor>
      <arglist>(TSDC &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ASGapplyAdjFO.html</anchorfile>
      <anchor>a21fd9a20aa38f6edf8e2f3d19d8d3c83</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASGapplyAdjTangentFO</name>
    <filename>classRVL_1_1ASGapplyAdjTangentFO.html</filename>
    <base>RVL::TSFO</base>
    <member kind="function">
      <type></type>
      <name>ASGapplyAdjTangentFO</name>
      <anchorfile>classRVL_1_1ASGapplyAdjTangentFO.html</anchorfile>
      <anchor>ae83ac42a3a5f6fd0acde6867b60a80e9</anchor>
      <arglist>(ASGaux const &amp;_aux, bool const &amp;_integerstep)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ASGapplyAdjTangentFO</name>
      <anchorfile>classRVL_1_1ASGapplyAdjTangentFO.html</anchorfile>
      <anchor>aa517a83a16b525259aa3f238d6d168ad</anchor>
      <arglist>(ASGapplyAdjTangentFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASGapplyAdjTangentFO</name>
      <anchorfile>classRVL_1_1ASGapplyAdjTangentFO.html</anchorfile>
      <anchor>abe37ac773dace266bbcdcd094735df30</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ASGapplyAdjTangentFO.html</anchorfile>
      <anchor>aa29c7ce4e691728f453320121138a981</anchor>
      <arglist>(TSDC &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ASGapplyAdjTangentFO.html</anchorfile>
      <anchor>a66574d72eba04036b6e5fa9a4b039850</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASGapplyFO</name>
    <filename>classRVL_1_1ASGapplyFO.html</filename>
    <base>RVL::TSFO</base>
    <member kind="function">
      <type></type>
      <name>ASGapplyFO</name>
      <anchorfile>classRVL_1_1ASGapplyFO.html</anchorfile>
      <anchor>a70a90b7049a9be32d93e256784516c81</anchor>
      <arglist>(ASGaux const &amp;_aux, bool const &amp;_integerstep)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ASGapplyFO</name>
      <anchorfile>classRVL_1_1ASGapplyFO.html</anchorfile>
      <anchor>ac78315c0d86f7595f6bead084c2e1684</anchor>
      <arglist>(ASGapplyFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASGapplyFO</name>
      <anchorfile>classRVL_1_1ASGapplyFO.html</anchorfile>
      <anchor>a5e23a12d9ba893ce5092583b848037da</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ASGapplyFO.html</anchorfile>
      <anchor>af853f8d036d2c2721dc859394d34fd5e</anchor>
      <arglist>(TSDC &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ASGapplyFO.html</anchorfile>
      <anchor>a2e0a1bc221b839d856638744d476ce78</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASGapplyTangentFO</name>
    <filename>classRVL_1_1ASGapplyTangentFO.html</filename>
    <base>RVL::TSFO</base>
    <member kind="function">
      <type></type>
      <name>ASGapplyTangentFO</name>
      <anchorfile>classRVL_1_1ASGapplyTangentFO.html</anchorfile>
      <anchor>a32de70550b7356cc6db78fda9bcb7c27</anchor>
      <arglist>(ASGaux const &amp;_aux, bool const &amp;_integerstep)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ASGapplyTangentFO</name>
      <anchorfile>classRVL_1_1ASGapplyTangentFO.html</anchorfile>
      <anchor>ae4577db1a01db92fd9ea17cd036f449e</anchor>
      <arglist>(ASGapplyTangentFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASGapplyTangentFO</name>
      <anchorfile>classRVL_1_1ASGapplyTangentFO.html</anchorfile>
      <anchor>ab828f7a2d9ad5e3d9bf9ab77aa83eb0f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ASGapplyTangentFO.html</anchorfile>
      <anchor>ae80aad621334f40dd22c18f188d63c31</anchor>
      <arglist>(TSDC &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ASGapplyTangentFO.html</anchorfile>
      <anchor>a14b10d1211a129fbd8906cd08ec4fb71</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASGaux</name>
    <filename>classRVL_1_1ASGaux.html</filename>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>ASGaux</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>aa3c124e2df9a88183c78b4134bbbcf03</anchor>
      <arglist>(GridDomain const &amp;gdom, int maxoff, float dt, std::vector&lt; int &gt; nls, std::vector&lt; int &gt; nrs, float amp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ASGaux</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>afcb22cac457514721415c76e3c2153e7</anchor>
      <arglist>(ASGaux const &amp;aux)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASGaux</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a2e656bc236a43620e61b11fcdbcbb0e9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>ae0e9dde6ecad51ea09f9646ac46539e6</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="variable">
      <type>GridDomain const &amp;</type>
      <name>gdom</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>afdbf8d41e785b2ebceea400dbbf6c988</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>dt</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>ad8db62c81bea775affc6603cff2b088e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>amp</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>af7c0edcaa30d354218c718565c93a4f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; int &gt;</type>
      <name>nls</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>aa19699c2af3afc247efd4d1291b6df8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; int &gt;</type>
      <name>nrs</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a06998d10b7912a3e7495c336329d62d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>lbc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>ae8b5cd2036edbd92614b2e860745a910</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>rbc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a232ee8d855e1ba10059fd28b9a634959</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float *</type>
      <name>vdiv_alloc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>ae0e54bc54e49c507f7a749e39301d529</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float *</type>
      <name>vdiv</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a3848c98760def7d31220e138c1af26f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>pgrad_alloc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a9298419df5808ccae580cac7b515d7a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>pgrad</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a693a09d8cb6892027944479d61134cd6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>ep</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>ad2e60903aec7220ecde3a7a9106f4959</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>epp</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a32d689c90a792d8288676328db040006</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>ev</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a3713f8a4551efbe3e5d71a11770a8489</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>evp</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a4491f7ba581bfcd1c812aa2fe9f27484</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>ep_alloc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>acef2af8662f87567e5d6adf8bffe6984</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>epp_alloc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a0a216bdf265dc7687a5342ed388d8f21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>ev_alloc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>abd8690c3660c5aadd5eb4a5892c9f6d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>evp_alloc</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>aeca23743a492207a3359de37c614c644</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float **</type>
      <name>coeffs</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>ae49ef92b538a6dec790c897122708c44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>maxoff</name>
      <anchorfile>classRVL_1_1ASGaux.html</anchorfile>
      <anchor>a2afc6a7c6dd93642031bb9e153a7064b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASGStep</name>
    <filename>classRVL_1_1ASGStep.html</filename>
    <base>TSStep&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>ASGStep</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>ad0afb196c85df479d81128864b455b74</anchor>
      <arglist>(ASGStep const &amp;stp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ASGStep</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>af16eacb3ff202bd6a11e36ce662c6f58</anchor>
      <arglist>(Grid const &amp;phys, int maxoff, float _tmin, float _dt, std::vector&lt; int &gt; _nls, std::vector&lt; int &gt; _nrs, float amp)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getTimeStep</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>ac47bb0136fed5640f3c561b7aa836b99</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; float &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>aea79d4be55147839de0ff56992c71b02</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; float &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>ac43a2905ab1f9648bcff88543f171895</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>ae16cb710c7518fab6ca70007a6b1c912</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>stepTimeFwd</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>a622c7da700edf5a279f7dcce8030f2f5</anchor>
      <arglist>(float &amp;t) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>stepTimeBwd</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>ae5bc10573d1c17b6f053f258439a8c4d</anchor>
      <arglist>(float &amp;t) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; float &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ASGStep.html</anchorfile>
      <anchor>a3cfcd5fc1e66e66dc89e47bbb711e496</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Axis</name>
    <filename>classRVL_1_1Axis.html</filename>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>Axis</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>aca2d7c1526c3e883baf32e1a215052dd</anchor>
      <arglist>(size_t _n=0, float _d=1.0f, float _o=0.0f, int _id=0, std::string _dunit=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Axis</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>ac649c88a9d4f03708f39c9c68f9d96e4</anchor>
      <arglist>(Axis const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a2d37649f2c7d5573941a663876c74e92</anchor>
      <arglist>(Axis const &amp;a) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator=</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a5785093074d97e63dca186613c014e1e</anchor>
      <arglist>(Axis const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLength</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a591edcabf6b9130515bae0bc71d081a9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getStep</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a46a963666e873af421ea97427bc1dcc3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getOrigin</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a2f8c08e41282ce5b1756161b195d3ccd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getID</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a720eb4ecfcd663fd436005f393ce5c35</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getUnit</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a207d8e1dc8de0e3fa104532451971e27</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getAxisType</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>afd345e7702a1198777f9960896033d7a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getAxisStartIndex</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a2a6df940119b6969cd8712eb427bc852</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getAxisEndIndex</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a4b288b875a7c54c36670f02c005b5b73</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getSubAxisStartIndex</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>ad10525576879051ab013a57037a61d18</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getSubAxisEndIndex</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a5c0dcca48dfcb77f6917320d1b1ab25d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSubAxisStartIndex</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>ac78a7360e8bb1a21c23a63d04e0ada23</anchor>
      <arglist>(int _gs) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSubAxisEndIndex</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>a4635ad847f725ddb34b83521bf4f84c5</anchor>
      <arglist>(int _ge) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getAxisSize</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>aa07c7a42b6411091008c57eefb0851aa</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Axis.html</anchorfile>
      <anchor>ad7ea4c2758bd1787224533c0e81c1e92</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Grid</name>
    <filename>classRVL_1_1Grid.html</filename>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>Grid</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a5eece1b45ee93a3c543b9ba24e34b6db</anchor>
      <arglist>(size_t _gdim=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Grid</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>af668045c808286ece38f27ff3af72dbc</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>addAxis</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>ac3954682d948613b01cf7a5975a5e550</anchor>
      <arglist>(Axis const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDimension</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a1088cfe1fbb3af778337861837f0b553</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isInit</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a9c157927411278365a8e360e6921ad13</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Axis const &amp;</type>
      <name>getAxis</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a6da62c04bf3d09b25387d012831cba21</anchor>
      <arglist>(size_t i) const </arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getCellVolume</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a73898e6c5bd5ee46135c49ae0ccd7340</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getAxisLimits</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>afd5d7e0d35a4cd950869f9f293a26543</anchor>
      <arglist>(IPNT gs0, IPNT ge0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getSubAxisLimits</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a07fb322c2c4cf831bf36f19de7a5ad31</anchor>
      <arglist>(IPNT gs, IPNT ge) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSubAxisLimits</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a58122b00d83092585a701d67b99cd81c</anchor>
      <arglist>(IPNT gs, IPNT ge) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a595962eec4fc70e519886510ff634051</anchor>
      <arglist>(Grid const &amp;g) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a1e2bb352e4fb81f1d5c63740cac146f4</anchor>
      <arglist>(Grid const &amp;g) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Grid.html</anchorfile>
      <anchor>a7e437a590a56d379e467f0010e60c268</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridCopyOverlapFO</name>
    <filename>classRVL_1_1GridCopyOverlapFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridCopyOverlapFO</name>
      <anchorfile>classRVL_1_1GridCopyOverlapFO.html</anchorfile>
      <anchor>a11e3da3fff188954f277d5bea36f8ec7</anchor>
      <arglist>(bool _plus=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridCopyOverlapFO</name>
      <anchorfile>classRVL_1_1GridCopyOverlapFO.html</anchorfile>
      <anchor>a134164444272550207b9b6530305f810</anchor>
      <arglist>(GridCopyOverlapFO const &amp;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridCopyOverlapFO</name>
      <anchorfile>classRVL_1_1GridCopyOverlapFO.html</anchorfile>
      <anchor>af40049bc2b51486650162670c97d313b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1GridCopyOverlapFO.html</anchorfile>
      <anchor>a65dadb8328fcdbff3e03de6aa9021027</anchor>
      <arglist>(LocalDataContainer&lt; float &gt; &amp;inside, LocalDataContainer&lt; float &gt; const &amp;outside)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1GridCopyOverlapFO.html</anchorfile>
      <anchor>a14dcbf36b452e019528921bcf4af375e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridDC</name>
    <filename>classRVL_1_1GridDC.html</filename>
    <member kind="function">
      <type></type>
      <name>GridDC</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>adb16d48825695df509268ecd6b94efb7</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridDC</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a30187d59fc32e5a5d5d3a9705ee6b63b</anchor>
      <arglist>(GridDC const &amp;d)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridDC</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>ad3085eb69b575bd16f94d0b44ba59201</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Grid const &amp;</type>
      <name>getGrid</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a64ebbf314dec55db9745c813c11a5895</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>getData1D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>afc7a502b4223b11169e77e4f909f3388</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const *</type>
      <name>getData1D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a7b948ec73242cb89ca86ac479b10daaf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>getData1DGlobal</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a2966a69390cd9feff922c59012d3b8b8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const *</type>
      <name>getData1DGlobal</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a03e400614935c549581aa354b8b84a64</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float **</type>
      <name>getData2D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>afdc2ea9bbae99a75d37d5c0074f97252</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const **</type>
      <name>getData2D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a73e3a3f1d8bf2cca4f41f641cb447af5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float **</type>
      <name>getData2DGlobal</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>ac02f4b37c888fa06c71e089d9f34db31</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const **</type>
      <name>getData2DGlobal</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a6b94349bb02ddd942012420798fe0c0b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float ***</type>
      <name>getData3D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>aab161088476ccd4a7d9673e505cd90cf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const ***</type>
      <name>getData3D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a41d584159bdea5bc9cb19a1d4fc0513f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float ***</type>
      <name>getData3DGlobal</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>abc621cacd10ce7f7276baace5cbc8941</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const ***</type>
      <name>getData3DGlobal</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a311daa63e070e963f2ff6620c5fc440b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float ****</type>
      <name>getData4D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a144ffc0cf8457c961044e298174a2a4e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const ****</type>
      <name>getData4D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>ae96eee327d82162d76827632f9a443f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>float *****</type>
      <name>getData5D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>a4bc77d26cbf0d41a9e1337a415cb4dcd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float const *****</type>
      <name>getData5D</name>
      <anchorfile>classRVL_1_1GridDC.html</anchorfile>
      <anchor>ab6f9f1ffc06a199542b4019a907a267d</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridDCFactory</name>
    <filename>classRVL_1_1GridDCFactory.html</filename>
    <base>RVL::DataContainerFactory</base>
    <member kind="function">
      <type></type>
      <name>GridDCFactory</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>ae9225591e7ca4abd4352783ea126be4c</anchor>
      <arglist>(Grid const &amp;_g)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridDCFactory</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>ada4d2a6bb283011a82f059f7af79eedd</anchor>
      <arglist>(GridDCFactory const &amp;gf)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>aaad1c4620e0db4d4c16f961dd467fcb2</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>afadb1a0f068ac0aee1edcd3aa2944be5</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>Grid const &amp;</type>
      <name>getGrid</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>af71eea7c86176748c13ea75c00c9f89f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>a9dc0977bf830f509aed5be5131e8efcc</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>GridDC *</type>
      <name>buildDC</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>a4534ccc0748304a5176fc1a2acb2e3cf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>DataContainer *</type>
      <name>build</name>
      <anchorfile>classRVL_1_1GridDCFactory.html</anchorfile>
      <anchor>ac194072cde09bfbdc596eba21fd8b421</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridDCIOFO</name>
    <filename>classRVL_1_1GridDCIOFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridDCIOFO</name>
      <anchorfile>classRVL_1_1GridDCIOFO.html</anchorfile>
      <anchor>a455c3532a4dc3d3452d57c8a7a270844</anchor>
      <arglist>(std::string _fname, bool _load=true, bool _extend=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridDCIOFO</name>
      <anchorfile>classRVL_1_1GridDCIOFO.html</anchorfile>
      <anchor>a9fd8958127591372bb0c2ddc013d4dd6</anchor>
      <arglist>(GridDCIOFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1GridDCIOFO.html</anchorfile>
      <anchor>a416b9a409f4cf45b6d1c04cf534b7d1f</anchor>
      <arglist>(LocalDataContainer&lt; float &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1GridDCIOFO.html</anchorfile>
      <anchor>a70791e371e45964fa083f89cfb34c5d5</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridDomain</name>
    <filename>classRVL_1_1GridDomain.html</filename>
    <base>TSSpace&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>GridDomain</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>a42916f2b6354f8fe5c266ded639c5493</anchor>
      <arglist>(std::vector&lt; Grid * &gt; ctrlgrids, std::vector&lt; Grid * &gt; stategrids, std::map&lt; std::string, int &gt; _ind)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridDomain</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>a29842ee1fd5f597193bfd67ae4975fd8</anchor>
      <arglist>(GridDomain const &amp;gd)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IPNT &gt; const &amp;</type>
      <name>get_gsa</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>aa7690e3a14ef976908f5796884bc8c61</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IPNT &gt; const &amp;</type>
      <name>get_gea</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>ad3e76ffb621c52a78b61909a91cfd32f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IPNT &gt; const &amp;</type>
      <name>get_gsc</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>abf560ab21b7a784897badc42897a904a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IPNT &gt; const &amp;</type>
      <name>get_gec</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>afffe61cd49c84a989ac2085e19ac13b1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; std::string, int &gt;</type>
      <name>get_ind</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>a9b7ebcbdc73d76f417a2a26c048ed6ae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getDimension</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>aa910b185dfed1cd5dff5479026ac14a0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; float &gt;</type>
      <name>getSteps</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>a3b3ba61067e73f5e65a6c3c26e81d7de</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1GridDomain.html</anchorfile>
      <anchor>a80b9599880021a25eed39f0df3096f4b</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridExtendFO</name>
    <filename>classRVL_1_1GridExtendFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridExtendFO</name>
      <anchorfile>classRVL_1_1GridExtendFO.html</anchorfile>
      <anchor>ac6975ba3be822e6a8d7671590ccac4a0</anchor>
      <arglist>(Grid const &amp;_g)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridExtendFO</name>
      <anchorfile>classRVL_1_1GridExtendFO.html</anchorfile>
      <anchor>aa49f70edbb89611fbec1d854816784c3</anchor>
      <arglist>(GridExtendFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1GridExtendFO.html</anchorfile>
      <anchor>a321c7ec501f080e1d5e81acaf72442c3</anchor>
      <arglist>(LocalDataContainer&lt; float &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1GridExtendFO.html</anchorfile>
      <anchor>a56c703a25177ae60775b183d04ede501</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridSpace</name>
    <filename>classRVL_1_1GridSpace.html</filename>
    <base>StdSpace&lt; float, float &gt;</base>
    <member kind="function">
      <type></type>
      <name>GridSpace</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>afce4297cae695ea4da86bf073c44fdde</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridSpace</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>a95e1555fb6451a0595a8da47e36650db</anchor>
      <arglist>(GridSpace const &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type>LinearAlgebraPackage&lt; float &gt; const &amp;</type>
      <name>getLAP</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>a453a2077392d0df1231b69bc4a8c40db</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataContainerFactory const &amp;</type>
      <name>getDCF</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>a41a3e62d018a6f57d863cb4a40d64f39</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Grid const &amp;</type>
      <name>getGrid</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>aff38cf7dc7478a006655b9b511041a5f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>a483749e2eeaf2c9d0854d4cce511a05c</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Space&lt; float &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1GridSpace.html</anchorfile>
      <anchor>a406b0c70836b9af629ee39f17af1f0a6</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GridtoTSOp</name>
    <filename>classRVL_1_1GridtoTSOp.html</filename>
    <base>TSSample&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>GridtoTSOp</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>aecd0c10e354bf71f0de4c17832fca265</anchor>
      <arglist>(float _tsamp, std::vector&lt; size_t &gt; _ind, GridSpace const &amp;_gsp, Space&lt; float &gt; const &amp;_tsp, bool _extend=false, float _tol=100 *numeric_limits&lt; float &gt;::epsilon())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridtoTSOp</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a05953aed78f4d0e9898f7f8adb119b53</anchor>
      <arglist>(GridtoTSOp const &amp;gop)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTestTime</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>af4e0a12783f83f8dd6993c9254acae9a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; float &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a8dc3601be9efb2882355304a4b2757ed</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; float &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a7de931ad74012b85afd3864d6bbf2292</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyPlusOp</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>ac5cf685d8b292eaa4770814bc4229711</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyPlusAdjOp</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a5beda3bbb4d04b08f1bee804b0a853c5</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>ab429c073c24bbbd3bfaa8a1202e1cc48</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a486a79817f9fa068099bb4e58b5aeec4</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyPlus</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>ad82f207778954235789cb6149ccc0fb4</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a28314a938e2192b81c12af362a404428</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjPlus</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>a1a003f448196d2697bfc358e92740240</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; float &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1GridtoTSOp.html</anchorfile>
      <anchor>aba943de16ac906e7fccb735d67419ba7</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinRestrictTSStep</name>
    <filename>classRVL_1_1LinRestrictTSStep.html</filename>
    <templarg></templarg>
    <base>LinearOp&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>LinRestrictTSStep</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>aea6dc034ab177be7c87c334a9724a7e4</anchor>
      <arglist>(TSStep&lt; T &gt; const &amp;_step, Vector&lt; T &gt; const &amp;x0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinRestrictTSStep</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>a9d2880294cc544300f91304646cd8ac5</anchor>
      <arglist>(LinRestrictTSStep&lt; T &gt; const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>aa81911c614511cd2bcae1388de77acb6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>af3388d95ac071344176de46594721e46</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>a39b3ac27d5a2d983081ac95d7c5ee223</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>a7f0205323343e665212dd6cc53533422</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>a0d3d871ea572bea2ec7642ed46e49611</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LinRestrictTSStep.html</anchorfile>
      <anchor>a671cfeb53e9ecd03f10d63cf36712105</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TanRestrictTSStep</name>
    <filename>classRVL_1_1TanRestrictTSStep.html</filename>
    <templarg></templarg>
    <base>LinearOp&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>TanRestrictTSStep</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>a7a5753d5be76ef84d134b698527c520d</anchor>
      <arglist>(TSStep&lt; T &gt; const &amp;_step, Vector&lt; T &gt; const &amp;x0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TanRestrictTSStep</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>a47703d3cbef5fdb6e73a0e37f0fc2895</anchor>
      <arglist>(TanRestrictTSStep&lt; T &gt; const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>a66da0c02e02e4a67b1168a465ef17e8f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>aaa0b2d8a9fe2d759094a9cf6bbfd5429</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>a21295a8759004bee6908164f94daf910</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>af1ed982c66335b9722f6c8898fcd6c18</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>a496bc876395e8c1c000041a4634ae24b</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1TanRestrictTSStep.html</anchorfile>
      <anchor>a26a3b712252835c80eac18b0a041b428</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TimeStepOp</name>
    <filename>classRVL_1_1TimeStepOp.html</filename>
    <templarg></templarg>
    <base>LinOpValOp&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>TimeStepOp</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>ab8326f4d507222885aa80bbc3f05d670</anchor>
      <arglist>(TSStep&lt; T &gt; &amp;_step, TSSample&lt; T &gt; &amp;_src, TSSample&lt; T &gt; &amp;_data, bool _testflag=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TimeStepOp</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>ac179ab33167928c609cf0d2c1a224599</anchor>
      <arglist>(TimeStepOp&lt; T &gt; const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~TimeStepOp</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a0dc8b89b02cc3fdb139f407c43b7e35f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ProductSpace&lt; float &gt; const &amp;</type>
      <name>getProductDomain</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>ae05c32abd06ff4e5578658169ef161e5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; float &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a20cc144b200161bbce7e38e5ff7d1dff</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>testSrcAdj</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a4fa38ff123f7eb34b778eab9f8b74fd0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>testDataAdj</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a9ab646b9b2cd89d7b7c15e0a9f3c6ded</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>testStepAdj0</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>adcfc684f69f282fdfc125716e6dfa123</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;ctrl)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>testStepAdj1</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>aebdc6da67ae6308a3ed492914462481b</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;ctrl)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>testAll</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a7bb625998017e832d06be3b83d3235e4</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x0)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>ab7801aa3b187d90a57ff17ed22858363</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>nowandthen</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a957f821d87dc21f5f078d92190eabc86</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;ctrlstate, Vector&lt; T &gt; const &amp;x1, T tstart, T tstop) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>save</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>af6e4782fd1f60f74e5827ca3ea8401d1</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, T t) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>load</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>ac783702ada9ee34898741ffaefc23c85</anchor>
      <arglist>(Vector&lt; T &gt; &amp;x, T t) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply0</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a9e4b96786e2dbc995ff292f482f93af9</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x0, Vector&lt; T &gt; const &amp;x1, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj0</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a43ff40afa70720833431cd2a6c0f4718</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x0, const Vector&lt; T &gt; &amp;x1, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyPartialDeriv0</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a0749a0cdcd7ae40bea2b0ec459707e1c</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x0, const Vector&lt; T &gt; &amp;x1, const Vector&lt; T &gt; &amp;dx0, Vector&lt; T &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjPartialDeriv0</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>af5bda51f99066437893ad2d0a195e740</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x0, const Vector&lt; T &gt; &amp;x1, const Vector&lt; T &gt; &amp;dy, Vector&lt; T &gt; &amp;dx0) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>OperatorProductDomain&lt; T &gt; *</type>
      <name>clonePD</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>ac788b81013f39f1c7e8281764d32f149</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TSStep&lt; T &gt; &amp;</type>
      <name>step</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a7a928dd00e5e3da745d51f5e0327c808</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>TSSample&lt; T &gt; &amp;</type>
      <name>src</name>
      <anchorfile>classRVL_1_1TimeStepOp.html</anchorfile>
      <anchor>a40343752082f06d4d22a079910783e37</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TimeStepOpAllCache</name>
    <filename>classRVL_1_1TimeStepOpAllCache.html</filename>
    <templarg>T</templarg>
    <base>RVL::TimeStepOp</base>
    <member kind="function">
      <type></type>
      <name>TimeStepOpAllCache</name>
      <anchorfile>classRVL_1_1TimeStepOpAllCache.html</anchorfile>
      <anchor>a9bd8667c87043ec4e30ecd2fd99c9b1f</anchor>
      <arglist>(TSStep&lt; T &gt; &amp;_step, TSSample&lt; T &gt; &amp;_src, TSSample&lt; T &gt; &amp;_data, bool _testflag=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TimeStepOpAllCache</name>
      <anchorfile>classRVL_1_1TimeStepOpAllCache.html</anchorfile>
      <anchor>a5c69f3c1319dee3b004a7c61e653bcb8</anchor>
      <arglist>(TimeStepOpAllCache&lt; T &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TimeStepOpAllCache.html</anchorfile>
      <anchor>a938faf4994964c06a0ff3e40a9122ca8</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>save</name>
      <anchorfile>classRVL_1_1TimeStepOpAllCache.html</anchorfile>
      <anchor>a2fe26fc043f9723b62b1518ae920a230</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, T t) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>load</name>
      <anchorfile>classRVL_1_1TimeStepOpAllCache.html</anchorfile>
      <anchor>ac6e53acc5a985544bc89523a1b897940</anchor>
      <arglist>(Vector&lt; T &gt; &amp;x, T t) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinOpValOp&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1TimeStepOpAllCache.html</anchorfile>
      <anchor>aed0b8be29d8616bcd21d5efb747d0bb7</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TSDC</name>
    <filename>classRVL_1_1TSDC.html</filename>
    <base>RVL::StdProductDataContainer</base>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1TSDC.html</anchorfile>
      <anchor>a7fc3fae80164fcb53426f9fd9a110e85</anchor>
      <arglist>(FunctionObject &amp;f, std::vector&lt; DataContainer const * &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TSDC.html</anchorfile>
      <anchor>a3b1db889e4b6269108bc4ffa5063f19f</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TSFO</name>
    <filename>classRVL_1_1TSFO.html</filename>
    <base>RVL::FunctionObject</base>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1TSFO.html</anchorfile>
      <anchor>af70ab525e02948ee9c70ac46625a13b5</anchor>
      <arglist>(TSDC &amp;y) const  =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TSSample</name>
    <filename>classRVL_1_1TSSample.html</filename>
    <templarg>T</templarg>
    <base>LinearOp&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>TSSample</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a883402116f34744ea009c77085081eb6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>getMinTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>af127d30c9e3eb86d72deb012fb704d86</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>getMaxTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a8bc9c4ad7b456385583133b40f24ceb8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>getTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a1e313debc64ad12fee6acdc2323cbf88</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMinTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a61b88c885c457a47735877be41feafdf</anchor>
      <arglist>(T _tmin) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>ac476c5d97b7627be079b13a6d76c7c97</anchor>
      <arglist>(T _tmax) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a1fcc721702bec9de415c352bc9c49efd</anchor>
      <arglist>(T _t) const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setTestTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>ade14c5a7a0c08380e2ec5bbbe37ce9e7</anchor>
      <arglist>()=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSSample&lt; float &gt;</name>
    <filename>classRVL_1_1TSSample.html</filename>
    <base>LinearOp&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>TSSample</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a883402116f34744ea009c77085081eb6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getMinTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>af127d30c9e3eb86d72deb012fb704d86</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getMaxTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a8bc9c4ad7b456385583133b40f24ceb8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>getTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a1e313debc64ad12fee6acdc2323cbf88</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMinTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a61b88c885c457a47735877be41feafdf</anchor>
      <arglist>(float _tmin) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>ac476c5d97b7627be079b13a6d76c7c97</anchor>
      <arglist>(float _tmax) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>a1fcc721702bec9de415c352bc9c49efd</anchor>
      <arglist>(float _t) const</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setTestTime</name>
      <anchorfile>classRVL_1_1TSSample.html</anchorfile>
      <anchor>ade14c5a7a0c08380e2ec5bbbe37ce9e7</anchor>
      <arglist>()=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TSSpace</name>
    <filename>classRVL_1_1TSSpace.html</filename>
    <templarg>T</templarg>
    <base>StdProductSpace&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>TSSpace</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a74b7f15fc74ba43597092c30d06acec8</anchor>
      <arglist>(size_t nfac)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TSSpace</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>aeaad816e73b7a6d19306b2920046c8e8</anchor>
      <arglist>(TSSpace&lt; T &gt; const &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a416684490aa3c43b87d104c9f8ddb95f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a0617e8400895a85354458be37e65d8cb</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Space&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a52b8127ef70e68f7948dc093311d336e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSSpace&lt; float &gt;</name>
    <filename>classRVL_1_1TSSpace.html</filename>
    <base>StdProductSpace&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>TSSpace</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a74b7f15fc74ba43597092c30d06acec8</anchor>
      <arglist>(size_t nfac)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TSSpace</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>aeaad816e73b7a6d19306b2920046c8e8</anchor>
      <arglist>(TSSpace&lt; float &gt; const &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a416684490aa3c43b87d104c9f8ddb95f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a0617e8400895a85354458be37e65d8cb</anchor>
      <arglist>(ostream &amp;str) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Space&lt; float &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1TSSpace.html</anchorfile>
      <anchor>a52b8127ef70e68f7948dc093311d336e</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TSStep</name>
    <filename>classRVL_1_1TSStep.html</filename>
    <templarg>T</templarg>
    <base>Operator&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>TSStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a69c350d551261dd49d58b9c0e91497af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TSStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>abae4956f94eeb95130a611d17ad1c5eb</anchor>
      <arglist>(TSStep&lt; T &gt; &amp;ts)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TSStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>aeaa85306b0a733560e4aaddaf42a61fb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>TSSpace&lt; T &gt; const &amp;</type>
      <name>getTSDomain</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a8a68b6aeb3edfcec12bb99e4a8320f12</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>af88159a03c8d63d54cb3a10b2d14914d</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a974890a11d713fbbe28b7d881888eba1</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyTangent</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a466d5bab2fb8a85580a24d8b8642111a</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;xdx, Vector&lt; T &gt; &amp;ydy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjTangent</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>abb8c061e147b632b4546e36fd8f1c5dd</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;xdx, Vector&lt; T &gt; &amp;ydy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>addeb6db7cc456b59d1b56578548a097c</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, Vector&lt; T &gt; const &amp;dx, Vector&lt; T &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a79a2a52ccf6118e7fdee897b8e9d58f7</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, Vector&lt; T &gt; const &amp;dx, Vector&lt; T &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual T</type>
      <name>getTimeStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a1b9e1b3f56ad1f18b0b0333a96906b07</anchor>
      <arglist>() const  =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>stepTimeFwd</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a37e9d45eae66b7f12b41f5f5265b98f9</anchor>
      <arglist>(T &amp;t) const  =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>stepTimeBwd</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a0d18fe89ecd28c16dcf7961cbe2e1172</anchor>
      <arglist>(T &amp;t) const  =0</arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>TimeStepOp&lt; T &gt;</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>ae93a6a905d97be5378babc2ebe31b83d</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>LinRestrictTSStep&lt; T &gt;</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a995651945383dd73f7b6c71afeecbd58</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>TanRestrictTSStep&lt; T &gt;</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a3ba647cd6c20b2cb39e5c1ac43ca1072</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSStep&lt; float &gt;</name>
    <filename>classRVL_1_1TSStep.html</filename>
    <base>Operator&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>TSStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a69c350d551261dd49d58b9c0e91497af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TSStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>abae4956f94eeb95130a611d17ad1c5eb</anchor>
      <arglist>(TSStep&lt; float &gt; &amp;ts)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TSStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>aeaa85306b0a733560e4aaddaf42a61fb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>TSSpace&lt; float &gt; const &amp;</type>
      <name>getTSDomain</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a8a68b6aeb3edfcec12bb99e4a8320f12</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>af88159a03c8d63d54cb3a10b2d14914d</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a974890a11d713fbbe28b7d881888eba1</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; &amp;y) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyTangent</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a466d5bab2fb8a85580a24d8b8642111a</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;xdx, Vector&lt; float &gt; &amp;ydy) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjTangent</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>abb8c061e147b632b4546e36fd8f1c5dd</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;xdx, Vector&lt; float &gt; &amp;ydy) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>addeb6db7cc456b59d1b56578548a097c</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; const &amp;dx, Vector&lt; float &gt; &amp;dy) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a79a2a52ccf6118e7fdee897b8e9d58f7</anchor>
      <arglist>(Vector&lt; float &gt; const &amp;x, Vector&lt; float &gt; const &amp;dx, Vector&lt; float &gt; &amp;dy) const</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual float</type>
      <name>getTimeStep</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a1b9e1b3f56ad1f18b0b0333a96906b07</anchor>
      <arglist>() const  =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>stepTimeFwd</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a37e9d45eae66b7f12b41f5f5265b98f9</anchor>
      <arglist>(float &amp;t) const  =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>stepTimeBwd</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a0d18fe89ecd28c16dcf7961cbe2e1172</anchor>
      <arglist>(float &amp;t) const  =0</arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>TimeStepOp&lt; T &gt;</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>ae93a6a905d97be5378babc2ebe31b83d</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>LinRestrictTSStep&lt; T &gt;</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a995651945383dd73f7b6c71afeecbd58</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>TanRestrictTSStep&lt; T &gt;</name>
      <anchorfile>classRVL_1_1TSStep.html</anchorfile>
      <anchor>a3ba647cd6c20b2cb39e5c1ac43ca1072</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>RVL</name>
    <filename>namespaceRVL.html</filename>
    <class kind="class">RVL::ASGapplyAdjFO</class>
    <class kind="class">RVL::ASGapplyAdjTangentFO</class>
    <class kind="class">RVL::ASGapplyFO</class>
    <class kind="class">RVL::ASGapplyTangentFO</class>
    <class kind="class">RVL::ASGaux</class>
    <class kind="class">RVL::ASGStep</class>
    <class kind="class">RVL::Axis</class>
    <class kind="class">RVL::Grid</class>
    <class kind="class">RVL::GridCopyOverlapFO</class>
    <class kind="class">RVL::GridDC</class>
    <class kind="class">RVL::GridDCFactory</class>
    <class kind="class">RVL::GridDCIOFO</class>
    <class kind="class">RVL::GridDomain</class>
    <class kind="class">RVL::GridExtendFO</class>
    <class kind="class">RVL::GridSpace</class>
    <class kind="class">RVL::GridtoTSOp</class>
    <class kind="class">RVL::LinRestrictTSStep</class>
    <class kind="class">RVL::TanRestrictTSStep</class>
    <class kind="class">RVL::TimeStepOp</class>
    <class kind="class">RVL::TimeStepOpAllCache</class>
    <class kind="class">RVL::TSDC</class>
    <class kind="class">RVL::TSFO</class>
    <class kind="class">RVL::TSSample</class>
    <class kind="class">RVL::TSSpace</class>
    <class kind="class">RVL::TSStep</class>
    <member kind="function">
      <type>bool</type>
      <name>AdjointTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a4f4c54c14907a1a59176a38230bdd215</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;op, FunctionObject &amp;randomize, ostream &amp;str, int tol=100)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>DerivTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a54fd266fc449300086fc5a142b16bcba</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op, Vector&lt; Scalar &gt; const &amp;y, Vector&lt; Scalar &gt; const &amp;p, ostream &amp;str, int n=10, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType hmin=0.1, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType hmax=1.0, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType minrat=1.95)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>GradientTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ab5912f017113d3e671381e333e38d743</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;p, ostream &amp;str, int n=11, Scalar hmin=0.1, Scalar hmax=1.0, Scalar minrat=1.95)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataSize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a4718240d9213bb6f2393a2556eabf38e</anchor>
      <arglist>(MetaType const &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMetaSize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a81def3aac48080aee95fdc7ca91ec4bc</anchor>
      <arglist>(MetaType const &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Scan</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a1c248aea14e2b1280af6008324915a4e</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;p, int n=11, Scalar hmin=-ScalarFieldTraits&lt; Scalar &gt;::One(), Scalar hmax=ScalarFieldTraits&lt; Scalar &gt;::One(), ostream &amp;str=cout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SpaceTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>afa36714a2429f8ae2d8c292241918015</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;sp, Vector&lt; Scalar &gt; const &amp;v, std::string msg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>testRealOnly</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ac9553ebb5a1dfe83a513ff453a70e9ca</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ProtectedDivision</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a5e7beced6b8741d3583b4cb941b658ce</anchor>
      <arglist>(real a, real b, real &amp;quot, real tol=ScalarFieldTraits&lt; real &gt;::AbsZero())</arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>sgcoeffs</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a1f19b1c5eaccb7b96b056e291d7cbdd4</anchor>
      <arglist>(int k)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Grid * &gt;</type>
      <name>make_asg_ctrllist</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a9587e9bc6f5a0a3d5cc90d487741ad9d</anchor>
      <arglist>(Grid const &amp;phys, int maxoff)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Grid * &gt;</type>
      <name>make_asg_statelist</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a78994ac742904536267943832863507a</anchor>
      <arglist>(Grid const &amp;phys, int maxoff)</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; std::string, int &gt;</type>
      <name>make_asg_indices</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a6a4ebe8273d42505edd4481f1974478b</anchor>
      <arglist>(int dim)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>asg_pmlaxis</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a6b6b057b9b6bef3bf0158af1238ce840</anchor>
      <arglist>(int n0, int nl, int nr, float amp, float dt, int gtype, float **ep, float **epp)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>areCompatible</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a253d4433eb023895798e937d483a50f8</anchor>
      <arglist>(Axis const &amp;a, Axis const &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>areCompatible</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ad1d907ab39aa73cf5c50aa7f450b7ad2</anchor>
      <arglist>(Grid const &amp;g, Grid const &amp;h)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>writeMeta&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a153d53575b3e91a4d34c9cf4a4782f36</anchor>
      <arglist>(Grid const &amp;g, ostream &amp;e)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataSize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa35397ec153fa0b4c587e9186c0a96f7</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMetaSize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a8f67c9e5d41d97dbaa6aeccde5c8ca28</anchor>
      <arglist>(Grid const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>serialize&lt; Axis &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aec2738851d7d7dded515e4a99db99cbb</anchor>
      <arglist>(Axis const &amp;a, size_t &amp;len)</arglist>
    </member>
    <member kind="function">
      <type>Axis *</type>
      <name>deserialize&lt; Axis &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa278468793367e9344c5e43af61f1f1c</anchor>
      <arglist>(char *cbuf, size_t len)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>serialize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ae07fa7da4d6b2bbbb58b6705d9d43c0b</anchor>
      <arglist>(Grid const &amp;g, size_t &amp;len)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>deserialize&lt; Grid &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ae0bf0ea161cae5d14063865df8a0a357</anchor>
      <arglist>(char *cbuf, size_t len)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Grid &gt;</type>
      <name>make_padded_grid</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a434e8bfa1dc2bb99e22acdf49ab3265f</anchor>
      <arglist>(Grid const &amp;phys, std::vector&lt; int &gt; nlsloc, std::vector&lt; int &gt; nrsloc)</arglist>
    </member>
    <member kind="function">
      <type>T **</type>
      <name>dimView</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aef62f8bbac064b68be274998dc02857d</anchor>
      <arglist>(T *x, size_t n, size_t chunklen)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; Grid &gt;</type>
      <name>GridfromFile</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a167096aab0f1e00c5fdb9272ed2f5727</anchor>
      <arglist>(std::string fname)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>The RVL Time Stepping Operator Class</title>
    <filename>index</filename>
  </compound>
</tagfile>
