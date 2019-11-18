my $expected_version = "1.";
print("import: Excel::Writer::XLSX\n");
use Excel::Writer::XLSX;

if (defined Excel::Writer::XLSX->VERSION) {
	my $given_version = Excel::Writer::XLSX->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart\n");
use Excel::Writer::XLSX::Chart;

if (defined Excel::Writer::XLSX::Chart->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Area\n");
use Excel::Writer::XLSX::Chart::Area;

if (defined Excel::Writer::XLSX::Chart::Area->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Area->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Area->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Bar\n");
use Excel::Writer::XLSX::Chart::Bar;

if (defined Excel::Writer::XLSX::Chart::Bar->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Bar->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Bar->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Column\n");
use Excel::Writer::XLSX::Chart::Column;

if (defined Excel::Writer::XLSX::Chart::Column->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Column->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Column->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Doughnut\n");
use Excel::Writer::XLSX::Chart::Doughnut;

if (defined Excel::Writer::XLSX::Chart::Doughnut->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Doughnut->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Doughnut->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Line\n");
use Excel::Writer::XLSX::Chart::Line;

if (defined Excel::Writer::XLSX::Chart::Line->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Line->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Line->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Pie\n");
use Excel::Writer::XLSX::Chart::Pie;

if (defined Excel::Writer::XLSX::Chart::Pie->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Pie->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Pie->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Radar\n");
use Excel::Writer::XLSX::Chart::Radar;

if (defined Excel::Writer::XLSX::Chart::Radar->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Radar->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Radar->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Scatter\n");
use Excel::Writer::XLSX::Chart::Scatter;

if (defined Excel::Writer::XLSX::Chart::Scatter->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Scatter->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Scatter->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chart::Stock\n");
use Excel::Writer::XLSX::Chart::Stock;

if (defined Excel::Writer::XLSX::Chart::Stock->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chart::Stock->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chart::Stock->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Chartsheet\n");
use Excel::Writer::XLSX::Chartsheet;

if (defined Excel::Writer::XLSX::Chartsheet->VERSION) {
	my $given_version = Excel::Writer::XLSX::Chartsheet->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Chartsheet->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Drawing\n");
use Excel::Writer::XLSX::Drawing;

if (defined Excel::Writer::XLSX::Drawing->VERSION) {
	my $given_version = Excel::Writer::XLSX::Drawing->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Drawing->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Examples\n");
use Excel::Writer::XLSX::Examples;

if (defined Excel::Writer::XLSX::Examples->VERSION) {
	my $given_version = Excel::Writer::XLSX::Examples->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Examples->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Format\n");
use Excel::Writer::XLSX::Format;

if (defined Excel::Writer::XLSX::Format->VERSION) {
	my $given_version = Excel::Writer::XLSX::Format->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Format->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::App\n");
use Excel::Writer::XLSX::Package::App;

if (defined Excel::Writer::XLSX::Package::App->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::App->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::App->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Comments\n");
use Excel::Writer::XLSX::Package::Comments;

if (defined Excel::Writer::XLSX::Package::Comments->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Comments->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Comments->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::ContentTypes\n");
use Excel::Writer::XLSX::Package::ContentTypes;

if (defined Excel::Writer::XLSX::Package::ContentTypes->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::ContentTypes->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::ContentTypes->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Core\n");
use Excel::Writer::XLSX::Package::Core;

if (defined Excel::Writer::XLSX::Package::Core->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Core->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Core->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Custom\n");
use Excel::Writer::XLSX::Package::Custom;

if (defined Excel::Writer::XLSX::Package::Custom->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Custom->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Custom->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Packager\n");
use Excel::Writer::XLSX::Package::Packager;

if (defined Excel::Writer::XLSX::Package::Packager->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Packager->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Packager->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Relationships\n");
use Excel::Writer::XLSX::Package::Relationships;

if (defined Excel::Writer::XLSX::Package::Relationships->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Relationships->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Relationships->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::SharedStrings\n");
use Excel::Writer::XLSX::Package::SharedStrings;

if (defined Excel::Writer::XLSX::Package::SharedStrings->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::SharedStrings->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::SharedStrings->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Styles\n");
use Excel::Writer::XLSX::Package::Styles;

if (defined Excel::Writer::XLSX::Package::Styles->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Styles->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Styles->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Table\n");
use Excel::Writer::XLSX::Package::Table;

if (defined Excel::Writer::XLSX::Package::Table->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Table->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Table->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::Theme\n");
use Excel::Writer::XLSX::Package::Theme;

if (defined Excel::Writer::XLSX::Package::Theme->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::Theme->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::Theme->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::VML\n");
use Excel::Writer::XLSX::Package::VML;

if (defined Excel::Writer::XLSX::Package::VML->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::VML->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::VML->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Package::XMLwriter\n");
use Excel::Writer::XLSX::Package::XMLwriter;

if (defined Excel::Writer::XLSX::Package::XMLwriter->VERSION) {
	my $given_version = Excel::Writer::XLSX::Package::XMLwriter->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Package::XMLwriter->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Shape\n");
use Excel::Writer::XLSX::Shape;

if (defined Excel::Writer::XLSX::Shape->VERSION) {
	my $given_version = Excel::Writer::XLSX::Shape->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Shape->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Utility\n");
use Excel::Writer::XLSX::Utility;

if (defined Excel::Writer::XLSX::Utility->VERSION) {
	my $given_version = Excel::Writer::XLSX::Utility->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Utility->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Workbook\n");
use Excel::Writer::XLSX::Workbook;

if (defined Excel::Writer::XLSX::Workbook->VERSION) {
	my $given_version = Excel::Writer::XLSX::Workbook->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Workbook->VERSION . '
');

}
print("import: Excel::Writer::XLSX::Worksheet\n");
use Excel::Writer::XLSX::Worksheet;

if (defined Excel::Writer::XLSX::Worksheet->VERSION) {
	my $given_version = Excel::Writer::XLSX::Worksheet->VERSION;
	$given_version =~ s/0+$//;
	die('Expected version ' . $expected_version . ' but found ' . $given_version) unless ($expected_version eq $given_version);
	print('	using version ' . Excel::Writer::XLSX::Worksheet->VERSION . '
');

}
