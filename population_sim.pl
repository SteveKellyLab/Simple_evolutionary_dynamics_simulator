#! /usr/bin/perl -w

use Math::Random::MT::Auto qw(binomial);

&run_program;

sub run_program	{				
				&define_variables;
				system "clear";	
				print "Running Evolutionary Dynamics Simulation\n";				
				print "Generations\t$generations\n";
				print "Replicates\t$replicates\n";
				print "Heterozygous\t$h\n";
				print "Selection\t$s\n";	
				print "Pop. sizes\t$print_pop\n";			
				foreach my $N (@pop_sizes)
					{					
					&run_simulation($N);
					}
				&analyse_data;
				&print_data;
				}

sub define_variables	{						
						$replicates = 100;
						$generations = 30000000;						
						$s = 0.01;
						$h = 0.5;
						$print_freq = 1000;
						@pop_sizes = (10000, 100000, 1000000);
						$print_pop = join "\, ", @pop_sizes;
						%fixation_time = ();
						my $t = 0;
						while ($t <= $generations)
							{
							$time{$t} = 1;
							$t += $print_freq;
							}
						@sort_t = sort {$a <=> $b} keys %time;
						}

sub run_simulation	{
					my $N = $_[0];
					my $i = 0;
					my $init_freq = 1/$N;
					while ($i < $replicates)
						{
						my $af = $init_freq;
						my $x = 0;
						my $n = 2*$N;
						my $g = 0;						
						while ($g < $generations)
							{		
							my $temp = ((($af**2)*(1+$s)) + (($af*(1-$af))*(1+$h*$s)))/(1 + (2*($af*(1-$af)*$h*$s)) + (($af**2)*$s));
							my $x = (binomial($temp, $n))/($n);
							if ($x == 0)
								{
								$x = $af;
								}
							$af = $x;		
							++$g;
							if ($af == 1)
								{
								$fixation_time{$N}{$i} = $g;
								last;
								}
							}
						++$i;
						}
					}

sub analyse_data	{
					foreach my $N (@pop_sizes)
						{							
						foreach my $t (@sort_t)
							{												
							if (keys %{$fixation_time{$N}} == 0)
								{
								$data{$t}{$N}{'fixp'} = 0;
								$data{$t}{$N}{'fixt'} = 'NA';
								$data{$t}{$N}{'fixt_sd'} = 'NA';								
								}
							else
								{								
								my %temp = ();
								my $i = 0;
								my @sort = sort {$a <=> $b} values %{$fixation_time{$N}};								
								foreach my $e (@sort)
									{
									if ($e <= $t)
										{
										$temp{$i} = $e;
										++$i;
										}
									else
										{
										last;
										}
									}
								
								if (keys %temp == 0)
									{
									$data{$t}{$N}{'fixp'} = 0;
									$data{$t}{$N}{'fixt'} = 'NA';
									$data{$t}{$N}{'fixt_sd'} = 'NA';
									}
								else
									{
									my $pfix = (keys %temp)/($replicates);
									my $fix_c = 0;
									my $fix_t = 0;
									my $fix_sd = 0;
									foreach my $e (keys(%temp))
										{
										$fix_t += $temp{$e};
										++$fix_c;
										}

									my $tfix = $fix_t/$fix_c;

									foreach my $e (keys(%temp))
										{
										$fix_sd += ($temp{$e} - $tfix)**2;
										++$fix_c;
										}

									$fix_sd = int(sqrt($fix_sd/($fix_c - 1)));
									$tfix = int($tfix);
									
									$data{$t}{$N}{'fixp'} = $pfix;
									$data{$t}{$N}{'fixt'} = $tfix;
									$data{$t}{$N}{'fixt_sd'} = $fix_sd;									
									}
								}
							}
						}
					}
					
sub print_data	{
				open OUT, ">Simulation_result_s_$s\_h\_$h\_g_$generations\_rep_$replicates\.txt\n";
				print OUT "Generation";
				foreach my $N (@pop_sizes)
					{
					print OUT "\t$N\_pfix\t$N\_tfix\t$N\_tfix_sd";
					}
				print OUT "\n";
				foreach my $t (@sort_t)
					{
					print OUT "$t";
					foreach my $N (@pop_sizes)
						{
						print OUT "\t$data{$t}{$N}{'fixp'}\t$data{$t}{$N}{'fixt'}\t$data{$t}{$N}{'fixt_sd'}";
						}
					print OUT "\n";
					}
				close OUT;				
				}

exit;










