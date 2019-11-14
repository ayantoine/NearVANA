#!/usr/bin/perl -w

# Authors: Linux de Folis / Denis Filloux - CIRAD - UMR BGPI - 2017

# http://search.cpan.org/~jmcnamara/Excel-Writer-XLSX-0.85/

# Usage: ./Tab2Xls.pl 1_Tabular_file.tab (...) N_Tabular_file.tab Output.xlsx Sheetname_lentgh
# Output: Output.xlsx

use strict;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;

print "./Tab2Xls.pl @ARGV\n";

my $workbook = Excel::Writer::XLSX->new($ARGV[@ARGV-2]);
my $format = $workbook->add_format();
$format->set_num_format('@');
$format->set_align('left');
$format->set_font('Arial');
$format->set_size('10');

my @i = ();
foreach my $k (0..@ARGV-3) {
	my $sheetname = substr($ARGV[$k], 0, $ARGV[@ARGV-1]);
	if (length($sheetname) > 31) {
		print "Be carefull : the length of the sheetname $sheetname exceeds 31 characters !\n";
		$sheetname = substr($sheetname, 0, 15)."(...)".substr($sheetname, length($sheetname)-11, 31);
	}
	push (@i, $sheetname);
}

my $j = 0;
my $col = 0;
foreach (@i) {
	my $worksheet = "worksheet$_";
	$worksheet = $workbook->add_worksheet($_);
	$worksheet->set_zoom(90);
	open (TABFILE, $ARGV[$j]) or die "$ARGV[$j] : $!";
	my $row = 0;
	while (<TABFILE>) {
		chomp;
		my @Fld = split('\t', $_);
		$col = 0;
		foreach my $token (@Fld) {
			$worksheet->write_string($row, $col, $token, $format);
			$col++;
			if (length($token) >= 32768)
			{
			my $nbcar = length($token);
			my $cell = xl_rowcol_to_cell($row,$col);
			print "Be carefull : the cell $ARGV[$j]: $cell exceeds 32767 characters: $nbcar !\n";
			}
		}
		$row++;
	}
	if ($row >= 1048577)
	{
		print "Be carefull : the file $ARGV[$j] exceeds 1048576 lines: $row !\n";
	}
	$worksheet->autofilter( 0, 0, 0, $col-1);
	$worksheet->freeze_panes('D2'); # amelioration: inclure dans les arguments quelle cellule a freezer  	
	close (TABFILE);
	$j++;
}
