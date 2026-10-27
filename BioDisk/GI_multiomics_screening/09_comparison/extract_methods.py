from pathlib import Path
import zipfile,xml.etree.ElementTree as E
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/A_CRLM_phosphoregulation')
for p in (R/'02_raw/A_supp').glob('*.docx'):
 with zipfile.ZipFile(p) as z:
  t=E.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
  paras=[''.join(y.text or '' for y in x.findall('.//w:t',ns)) for x in t.findall('.//w:p',ns)]
  (R/'00_docs/supp_methods.txt').write_text('\n'.join(paras))
  import re
  for s in paras:
   if re.search('normaliz|imput|localiz|batch|TMT|ratio|MaxQuant',s,re.I):print(s)
