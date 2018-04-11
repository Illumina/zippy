# -*- mode: python -*-

#just an example.  With full paths, this file can be used to create a binary using pyinstaller.
block_cipher = None


a = Analysis(['make_params.py'],
             pathex=['path/to/zippy'],
             binaries=[],
             datas=[('/path/to/python/lib/python2.7/site-packages/pyflow/pyflow.py', './pyflow'),
             ('/path/to/python/lib/python2.7/site-packages/pyflow/pyflow.py', './site-packages/pyflow')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

b = Analysis(['zippy.py'],
             pathex=['path/to/zippy'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
MERGE((a, 'make_params', 'make_params'), (b, 'zippy', 'zippy'))
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
pyz2 = PYZ(b.pure, b.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='make_params',
          debug=False,
          strip=False,
          upx=True,
          console=True )
exe2 = EXE(pyz2,
          b.scripts,
          exclude_binaries=True,
          name='zippy',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='make_params')
coll2 = COLLECT(exe2,
               b.binaries,
               b.zipfiles,
               b.datas,
               strip=False,
               upx=True,
               name='zippy')


