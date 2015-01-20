import unittest
from mock import patch,call,mock_open,Mock
from topolenum import topolenum as te


@patch('tempfile.NamedTemporaryFile')
class TestFIFOfileBaseClassTMPFILEclass(unittest.TestCase):
  
  def setUp(self):
    reload(te)
    self.TMPFILE = te.FIFOfile.TMPFILE
    self.TMPFILE.init_class(mode='b',wbuf=0,rbuf=0,suffix='',delete=True,
                            dir='/dummy/path',check_freq=5)
  
  def test_instance_counting_and_file_creation_and_naming(self,patched_NTF):
    
    tmpfile_instance1 = self.TMPFILE()
    self.assertEqual(tmpfile_instance1.instcount,self.TMPFILE.instcount)
    self.assertEqual(self.TMPFILE.instcount,1)
    patched_NTF.assert_called_with('rb',0,suffix='',prefix='FIFOfile001_',
                                   dir='/dummy/path',delete=True)
    self.TMPFILE()
    self.assertEqual(self.TMPFILE.instcount,2)
    patched_NTF.assert_called_with('rb',0,suffix='',prefix='FIFOfile002_',
                                   dir='/dummy/path',delete=True)
  
  def test_spooling_and_file_spool_interaction(self,patched_NTF):
    self.assertFalse(hasattr(self.TMPFILE,'file_spool'))
    self.TMPFILE.start_spooling()
    self.assertTrue(hasattr(self.TMPFILE,'file_spool'))
    self.assertEqual(len(self.TMPFILE.file_spool),0)
    tmpfile_instance = self.TMPFILE.spool()
    self.assertEqual(len(self.TMPFILE.file_spool),1)
    self.assertIs(tmpfile_instance,self.TMPFILE.pop_from_spool())
    self.assertEqual(len(self.TMPFILE.file_spool),0)
    self.assertIs(None,self.TMPFILE.pop_from_spool())
  
  def test_writing_handle(self,patched_NTF):
    tmpfile_instance = self.TMPFILE()
    self.assertFalse(hasattr(tmpfile_instance,'wh'))
    with patch('__builtin__.open',mock_open(),create=True) as patched_open:
      tmpfile_instance.open()
      patched_open.assert_called_once_with(patched_NTF.return_value.name,'wb',0)
      self.assertTrue(hasattr(tmpfile_instance,'wh'))
      self.assertIs(tmpfile_instance.wh,patched_open.return_value)
    tmpfile_instance.close()
    tmpfile_instance.wh.close.assert_called_once_with()
  
  @patch('os.path.getsize',side_effect=[10.0,20.0])
  def test_size_checking(self,patched_getsize,patched_NTF):
    tmpfile_instance = self.TMPFILE()
    self.assertSequenceEqual([tmpfile_instance.size for i in xrange(5)],
                             [0.0,0.0,0.0,0.0,0.0])
    patched_getsize.assert_not_called()
    self.assertEqual(tmpfile_instance.size,10.0)
    patched_getsize.assert_called_once_with(patched_NTF.return_value.name)
    self.assertSequenceEqual([tmpfile_instance.size for i in xrange(5)],
                             [10.0,10.0,10.0,10.0,10.0])
    patched_getsize.assert_called_once_with(patched_NTF.return_value.name)
    self.assertEqual(tmpfile_instance.size,20.0)
    self.assertSequenceEqual(patched_getsize.call_args_list,
                             [call(patched_NTF.return_value.name),
                              call(patched_NTF.return_value.name)],list)
  
  @patch('os.path.exists',return_value=False)
  def test_file_discarding(self,patched_exists,patched_NTF):
    tmpfile_instance = self.TMPFILE()
    self.assertIs(tmpfile_instance.rh,patched_NTF.return_value)
    tmpfile_instance.discard()
    tmpfile_instance.rh.close.assert_called_once_with()


@patch('os.path.exists',return_value=False)
@patch('os.path.realpath',return_value='/dummy/path')
@patch('topolenum.topolenum.TemporaryDirectory')
@patch('tempfile.NamedTemporaryFile',**{'return_value.name':'dummy_temp_file'})
class TestFIFOfileBaseClass(unittest.TestCase):
  
  def setUp(self):
    reload(te)
  
  @patch('topolenum.topolenum.FIFOfile.TMPFILE.init_class')
  def test_tmpdir_creation_and_TMPFILE_class_init(self,patched_TMPFILE_initcls,
                                                  patched_NTF,patched_TmpDir,
                                                  patched_realpath,
                                                  patched_exists):
    fifo_obj = te.FIFOfile(top_path='/dummy/top/path/arg',suffix='dummy_suffix')
    patched_TmpDir.assert_called_once_with(dir='/dummy/top/path/arg',
                                           prefix='FIFOworkspace_',
                                           suffix='dummy_suffix')
    patched_TMPFILE_initcls.assert_called_once_with('b',0,0,'dummy_suffix',True,
                                                    '/dummy/path',100)
  
  def test_starting_fifo_OUT_end(self,patched_NTF,patched_TmpDir,
                                 patched_realpath,patched_exists):
    fifo_obj = te.FIFOfile()
    self.assertFalse(hasattr(te.FIFOfile.TMPFILE,'file_spool'))
    patched_NTF.assert_not_called()
    self.assertFalse(hasattr(fifo_obj,'current_reading_file'))
    
    fifo_obj.start_OUT_end()
    patched_NTF.assert_called_once_with('rb',0,suffix='',prefix='FIFOfile001_',
                                        dir='/dummy/path',delete=True)
    self.assertSequenceEqual(te.FIFOfile.TMPFILE.file_spool,[],seq_type=list)
    self.assertEqual(fifo_obj.current_reading_file.name,'dummy_temp_file')
  
  def test_starting_fifo_IN_end(self,patched_NTF,patched_TmpDir,
                                patched_realpath,patched_exists):
    fifo_obj = te.FIFOfile()
    self.assertFalse(hasattr(fifo_obj,'current_writing_file'))
    
    fifo_obj.start_OUT_end()
    with patch('__builtin__.open',mock_open(),create=True) as patched_open:
      fifo_obj.start_IN_end()
      self.assertIs(fifo_obj.current_writing_file,fifo_obj.current_reading_file)
      patched_open.assert_called_once_with('dummy_temp_file','wb',0)
      self.assertIs(fifo_obj.current_writing_file.wh,patched_open.return_value)
  
  @patch('os.path.getsize',side_effect=[500.0,1100.0])
  @patch('__builtin__.open')#,new=mock_open(),create=True)
  def test_wh_retrieval(self,patched_open,patched_getsize,patched_NTF,
                        patched_TmpDir,patched_realpath,patched_exists):
    def patched_NTF_side_effect(*args,**kwargs):
      return_val = Mock()
      return_val.name = kwargs['prefix']
      return return_val
    patched_NTF.side_effect = patched_NTF_side_effect
    
    def assert_no_rollover():
      self.assertFalse(fifo_obj.current_reading_file.wh.close.called)
      self.assertFalse(patched_NTF.called)
      self.assertFalse(patched_open.called)
      self.assertIs(fifo_obj.current_writing_file,
                    fifo_obj.current_reading_file)
      self.assertSequenceEqual(fifo_obj.TMPFILE.file_spool,[],seq_type=list)
    
    # Start up FIFO and clear out calls to mocks
    fifo_obj = te.FIFOfile(max_file_size_GB=1000.0,size_check_freq=2)
    fifo_obj.start_OUT_end()
    fifo_obj.start_IN_end()
    patched_NTF.reset_mock()
    patched_open.reset_mock()
    
    # First and second calls should trigger nothing
    fifo_obj.wh
    fifo_obj.wh
    self.assertFalse(patched_getsize.called)
    assert_no_rollover()
    
    # Third call should trigger file size check ...
    fifo_obj.wh
    patched_getsize.assert_called_once_with('FIFOfile001_')
    patched_getsize.reset_mock()
    self.assertEqual(fifo_obj.current_writing_file._size,500.0)
    # ... but the file size should be insufficient to trigger a rollover
    assert_no_rollover()
    
    # Fourth and fifth calls should, again, trigger nothing
    fifo_obj.wh
    fifo_obj.wh
    self.assertFalse(patched_getsize.called)
    assert_no_rollover()
    
    # Finally, sixth call should trigger a second size check ...
    fifo_obj.wh
    patched_getsize.assert_called_once_with('FIFOfile001_')
    # ... and the returned file size should trigger a rollover to a new file, ...
    fifo_obj.current_reading_file.wh.close.assert_called_once_with()
    patched_NTF.assert_called_once_with('rb',0,prefix='FIFOfile002_',suffix='',
                                        dir='/dummy/path',delete=True)
    patched_open.assert_called_once_with('FIFOfile002_','wb',0)
    self.assertSequenceEqual(fifo_obj.TMPFILE.file_spool,
                             [fifo_obj.current_writing_file])
    # ... meaning the reading and writing files are now different ...
    self.assertIsNot(fifo_obj.current_writing_file,
                     fifo_obj.current_reading_file)
    self.assertEqual(fifo_obj.current_writing_file._size,0.0)
    self.assertEqual(fifo_obj.current_reading_file._size,1100.0)
    


if __name__ == "__main__":
  #import sys;sys.argv = ['', 'Test.testName']
  unittest.main()
