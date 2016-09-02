# -*- mode: python -*-
a = Analysis(['tigger.py'],
             pathex=['/Users/brett/Dropbox/Code/tigger'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='tigger',
          debug=False,
          strip=None,
          upx=True,
          console=True )
