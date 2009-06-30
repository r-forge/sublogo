<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
   "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head><title>Sublogo dendrograms</title></head>

<body>

<?php 
   $params=$_POST;
   if($params['seqs']){
     if($params['subsmat']){
       $job=time();
       $jobdir='jobs/'.$job.'/';
       mkdir($jobdir);

       $prefix=$jobdir . $job.'-seqs.';
       $seqfile=$prefix.'txt';
       $fh=fopen($seqfile,'w');
       fwrite($fh,$params['seqs']);
       fclose($fh);
      
       $cmd="PATH=/opt/src/users/tdhock/R-2.7.2/bin:/opt/csw/bin:/opt/ocf/bin:/usr/bin time ./argv.plot ".$seqfile.
       " ".escapeshellarg($params['subsmat']).
       " ".escapeshellarg($params['tit']).
       " ".escapeshellarg($params['subtit']).
       " ".escapeshellarg($params['cutline']).
       " ".escapeshellarg($params['dendwidth']).
       " ".escapeshellarg($params['plotwidth']).
       " ".escapeshellarg($params['plotheight']).
       " > ".$prefix."stdouterr 2>&1 ; ./sublogo-dendrogram-cleanup";
       $cmdfile=$prefix.'command';
       $fh=fopen($cmdfile,'w');
       fwrite($fh,$cmd);
       fclose($fh);
       exec($cmd);
       ?>
       <p><a href="<?php echo $jobdir ?>">
       <img alt="sublogo dendrogram" src="<?php echo $seqfile.'.png' ?>" />
       </a></p>
       <?php 
     }else{
       $subsmat_error="enter the name of a common substitution matrix.";
     }
   }else{
     $seqs_error="enter prealigned sequences, in fasta format.";
   }
?>

<?php
if(! $params['seqs']){
  $params = array('cutline' => 32,
		  'tit' => "Zinc finger protein recognition helix sequences",
		  'subtit' => "Selected to bind triplet GGC",
		  'plotwidth' => 6,
		  'plotheight' => 6,
		  'dendwidth' => 33,
		  'seqs' => ">DQGHRTR\nDQGHRTR\n>DVGHRSR\nDVGHRSR\n>ESGHLRR\nESGHLRR\n>ESSKRKR\nESSKRKR\n>SRRNLTR\nSRRNLTR\n>TKGYLYK\nTKGYLYK\n>PSGYLYK\nPSGYLYK\n>WTSRLKH\nWTSRLKH\n>DKGHLRR\nDKGHLRR\n>DGSHLKR\nDGSHLKR\n>DRSNLRK\nDRSNLRK\n>ERSKLTR\nERSKLTR\n>ERSKLSR\nERSKLSR\n",
		  'subsmat' => 'BLOSUM62');
}
?>


<form action="" method="post">

<input type="submit" value="make sublogo dendrogram" />

<table>
  <tr>
    <td>cutline</td>
    <td>
      <input 
	 type="text" 
	 name="cutline" 
	 value="<?php echo $params['cutline'] ?>"
	 />
    </td>
  </tr>
  <tr>
    <td>title</td>
    <td>
      <input 
	 type="text" 
	 name="tit" 
	 value="<?php echo $params['tit'] ?>"
	 />
    </td>
  </tr>
  <tr>
    <td>subtitle</td>
    <td>
      <input 
	 type="text" 
	 name="subtit" 
	 value="<?php echo $params['subtit'] ?>"
	 />
    </td>
  </tr>
  <tr>
    <td>inches tall</td>
    <td>
      <input 
	 type="text" 
	 name="plotheight" 
	 value="<?php echo $params['plotheight'] ?>"
	 />
    </td>
  </tr>
  <tr>
    <td>inches wide</td>
    <td>
      <input 
	 type="text" 
	 name="plotwidth" 
	 value="<?php echo $params['plotwidth'] ?>"
	 />
    </td>
  </tr>
  <tr>
    <td>dendrogram width %</td>
    <td>
      <input 
	 type="text" 
	 name="dendwidth" 
	 value="<?php echo $params['dendwidth'] ?>"
	 />
    </td>
  </tr>
  <tr>
    <td>sequences</td>
    <td><textarea 
	   name="seqs" 
	   rows="5" 
	   cols="15"
	   ><?php echo $params['seqs'] ?></textarea></td>
    <td><?php echo $seqs_error ?></td>
  </tr>
  <tr>
    <td>substitution matrix</td>
    <td><textarea 
	   name="subsmat" 
	   rows="5" 
	   cols="15"
	   ><?php echo $params['subsmat'] ?></textarea></td>
    <td><?php echo $subsmat_error ?></td>
  </tr>
</table>

</form>

<center>
<table>
<tr>
<td align="center">
    <a href="http://validator.w3.org/check?uri=referer"><img
        src="http://www.w3.org/Icons/valid-xhtml10"
        alt="Valid XHTML 1.0 Transitional" height="31" width="88" /></a>
</td>
</tr>
</table>
</center>

</body>

</html>
