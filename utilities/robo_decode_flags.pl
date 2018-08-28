#! /usr/bin/perl -w

my $in_flag = hex(shift);

my $header = "robo.h";

my $open = 0;
if ((exists($ENV{ROBO_SRC_DIR}))&&($open == 0)) {
    $open = 1;
    open(H,"$ENV{ROBO_SRC_DIR}/$header") || ($open = 0);
}
if ($open == 0) {
    $open = 1;
    open(H,"./$header") || ($open = 0);
}
if ($open == 0) {
    $open = 1;
    open(H,"../$header") || ($open = 0);
}
if ($open == 0) {
    $open = 1;
    open(H,"$ENV{HOME}/usr/src/astronomy/spectra/src/$header") || ($open = 0);
}

if ($open == 0) {
    die "Could not find robo.h header file.";
}

my %bitmasks;

my $scan = 0;
while (<H>) {
    chomp;
    if ($_ =~ /^\s*$/) {
	next;
    }
    if ($_ =~ /END LINE FLAG BITMASKS/) {
	$scan = 0;
    }
    if ($scan == 1) {
	my (undef,$name,$flag,@comment) = split /\s+/;
	my $comment = join ' ', @comment;
	$comment =~ s%/\* %%;
	$comment =~ s% \*/%%;
	$bitmasks{$name}{FLAG} = hex($flag);
	$bitmasks{$name}{COMMENT} = $comment;
    }
    if ($_ =~ /BEGIN LINE FLAG BITMASKS/) {
	$scan = 1;
    }
}
close(H);

foreach my $n (sort {$bitmasks{$a}{FLAG} <=> $bitmasks{$b}{FLAG}} (keys %bitmasks)) {
    if ($in_flag & $bitmasks{$n}{FLAG}) {
	printf("%-20s   0x%06x  %s\n",
	       $n,$bitmasks{$n}{FLAG},$bitmasks{$n}{COMMENT});
    }
}
	
